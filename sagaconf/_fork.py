"""Run independent report steps in forked children, with the output of a sequential run.

run_steps(steps, savedir, ...) does what `for step in steps: step(savedir)` does. Each step
runs in a forked child (Linux only) that writes into its own directory under savedir and,
instead of printing, records every write to sys.stdout/sys.stderr, every flush, every warning
shown and every change of the warning filters, in order. The steps are committed in step order:
the parent replays each child's record (warnings are issued again, so the parent's filters and
"once" registries decide what is shown, as in a sequential run) and moves its files into place.

A step that is not forked (no memory budget or CPU left, a failed fork) runs in-process in its
turn. A step whose child fails, or cannot be replayed exactly, runs again in-process after the
global state of the steps before it is restored, and so do all the steps after it: from the
first failure on, the run is the sequential one.

Memory: every child reserves an estimate of its peak extra memory from one budget shared by the
whole process tree (FRACTION of the cgroup or machine memory, less what other processes use);
a step that does not fit waits for a running child, or runs in-process if none is running.
"""
import contextlib
import fcntl
import io
import mmap
import os
import pickle
import shutil
import signal
import struct
import sys
import time
import warnings

FRACTION = 0.65  # of the memory limit (cgroup, else machine) that the process tree may use
SAVER_BYTES = 512 << 20  # a forked figure save
MAX_SAVERS = 4  # figure saves pending at once in one step's child

_S = {"pool": None, "live": {}, "dirs": set(), "parts": set(), "savers": [], "step_child": False,
      "depth": 0, "sig": None, "seq": 0, "rec": None}


def _cpus():
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


def usable():
    """True where steps can run in forked children (Linux, 2+ CPUs, standard warning display,
    library deprecation warnings ignored)."""
    return (sys.platform.startswith("linux") and hasattr(os, "fork") and hasattr(os, "memfd_create")
            and _cpus() >= 2 and all(hasattr(warnings, n) for n in (
                "_showwarnmsg", "_showwarnmsg_impl", "_showwarning_orig", "_filters_mutated"))
            and warnings.showwarning is warnings._showwarning_orig
            and getattr(warnings._showwarnmsg_impl, "__module__", None) == "warnings"
            and _deprecations_hidden())


def frame_bytes(*frames):
    """The bytes of the frames' columns (8 per object cell), from shapes and dtypes only."""
    return sum(len(df) * sum(dt.itemsize for dt in df.dtypes) for df in frames)


# ---- memory budget shared by the process tree ------------------------------------------------

def _read(path):
    try:
        with open(path) as fh:
            return fh.read()
    except OSError:
        return None


def _kb_field(text, key):
    for line in (text or "").splitlines():
        if line.startswith(key + ":"):
            return int(line.split()[1]) * 1024
    return None


def _usage(pid="self"):
    """Anonymous and shared memory resident in process pid."""
    status = _read("/proc/%s/status" % pid)
    anon = _kb_field(status, "RssAnon")
    if anon is None:
        anon = _kb_field(status, "VmRSS") or 0
    return anon + (_kb_field(status, "RssShmem") or 0)


def _cgroup():
    """(limit, usage) of the tightest memory cgroup this process is in, or (None, None)."""
    best = (None, None)
    for line in (_read("/proc/self/cgroup") or "").splitlines():
        hid, ctrl, path = line.split(":", 2)
        if hid == "0":  # cgroup v2
            parts = [p for p in path.split("/") if p]
            for k in range(len(parts), -1, -1):
                d = "/sys/fs/cgroup/" + "/".join(parts[:k])
                m = (_read(d + "/memory.max") or "").strip()
                if m.isdigit() and (best[0] is None or int(m) < best[0]):
                    stat = (_read(d + "/memory.stat") or "").splitlines()
                    best = (int(m), sum(int(l.split()[1]) for l in stat if l.split()[0] in ("anon", "shmem")))
        elif "memory" in ctrl.split(","):  # cgroup v1
            d = "/sys/fs/cgroup/memory" + path
            m = (_read(d + "/memory.limit_in_bytes") or "").strip()
            if m.isdigit() and int(m) < 1 << 60 and (best[0] is None or int(m) < best[0]):
                stat = (_read(d + "/memory.stat") or "").splitlines()
                best = (int(m), sum(int(l.split()[1]) for l in stat
                                    if l.split()[0] in ("total_rss", "total_shmem")))
    return best


def _budget():
    """Bytes the process tree may use: FRACTION of the memory limit, less what other processes
    of the same cgroup use now (this process's own use is subtracted at each reservation)."""
    own = _usage()
    avail = _kb_field(_read("/proc/meminfo"), "MemAvailable")
    total = None if avail is None else int(FRACTION * (avail + own))
    limit, usage = _cgroup()
    if limit is None:
        mb = os.environ.get("SLURM_MEM_PER_NODE")
        if not mb and os.environ.get("SLURM_MEM_PER_CPU"):
            mb = int(os.environ["SLURM_MEM_PER_CPU"]) * int(os.environ.get("SLURM_CPUS_ON_NODE", "1"))
        if mb:
            limit, usage = int(mb) << 20, own
    if limit is not None:
        cg = int(FRACTION * limit) - max(0, usage - own)
        total = cg if total is None else min(total, cg)
    return max(0, total or 0)


class _Pool:
    """Reservations (pid, bytes) in a shared memory file locked with lockf. A slot belongs to a
    forked child (pid > 0) or, while it forks, to its parent (-pid); slots of dead processes are
    freed when seen. The root process's own memory is counted as it is, not reserved."""
    SLOTS = 256

    def __init__(self):
        self.fd = os.memfd_create("sagaconf_budget")
        os.ftruncate(self.fd, 16 * self.SLOTS)
        self.mm = mmap.mmap(self.fd, 16 * self.SLOTS)
        self.root, self.total, self.cpus = os.getpid(), _budget(), _cpus()

    def reserve(self, nbytes):
        fcntl.lockf(self.fd, fcntl.LOCK_EX)
        try:
            used, n, free = 0, 0, None
            for i in range(self.SLOTS):
                pid, b = struct.unpack_from("qq", self.mm, 16 * i)
                if pid and not _alive(abs(pid)):
                    struct.pack_into("qq", self.mm, 16 * i, 0, 0)
                    pid = 0
                if not pid:
                    free = i if free is None else free
                    continue
                used, n = used + b, n + 1
            # the root process works too: at most cpus - 1 children
            if free is None or n + 1 >= self.cpus or used + nbytes + _usage(self.root) > self.total:
                return None
            struct.pack_into("qq", self.mm, 16 * free, -os.getpid(), nbytes)
            return free
        finally:
            fcntl.lockf(self.fd, fcntl.LOCK_UN)

    def assign(self, slot, pid):
        if slot is None:
            return
        fcntl.lockf(self.fd, fcntl.LOCK_EX)
        try:
            struct.pack_into("q", self.mm, 16 * slot, pid)
        finally:
            fcntl.lockf(self.fd, fcntl.LOCK_UN)


def _alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except OSError:
        pass
    return True


def _pool():
    if _S["pool"] is None:
        _S["pool"] = _Pool()
    return _S["pool"]


# ---- recording output ------------------------------------------------------------------------

class _Stream(io.TextIOBase):
    """Stands in for sys.stdout / sys.stderr in a child: records writes and flushes in order."""

    def __init__(self, rec, tag, real):
        self._rec, self._tag, self._real = rec, tag, real
        self._enc = getattr(real, "encoding", None)
        self._err = getattr(real, "errors", None) or "strict"

    @property
    def encoding(self):
        return self._enc

    @property
    def errors(self):
        return self._err

    def write(self, s):
        if not isinstance(s, str):
            raise TypeError("write() argument must be str, not %s" % type(s).__name__)
        if self._enc:
            s.encode(self._enc, self._err)  # fails where the real stream would
        if not self._rec.in_warning:
            self._rec.events.append((self._tag, s))
        return len(s)

    def flush(self):
        if not self._rec.in_warning:
            self._rec.events.append((self._tag + "f",))

    def writable(self):
        return True

    def isatty(self):
        return self._real.isatty()

    def fileno(self):
        return self._real.fileno()  # a write to it reaches the child's raw-output file

    def close(self):
        pass  # io.IOBase closes (and so flushes) an object it finalizes: a flush is output here

    def __del__(self):
        pass


class _Tee(_Stream):
    """Stands in for sys.stdout / sys.stderr in-process: records and passes everything on."""

    def write(self, s):
        n = self._real.write(s)
        if not self._rec.in_warning:
            self._rec.events.append((self._tag, s))
        return n

    def flush(self):
        self._real.flush()
        if not self._rec.in_warning:
            self._rec.events.append((self._tag + "f",))

    def fileno(self):
        self._rec.exact = False  # output written to the descriptor is not recorded
        return self._real.fileno()

    def __getattr__(self, name):
        self._rec.exact = False
        return getattr(self._real, name)


class _Recorder:
    """Records writes to sys.stdout/sys.stderr, flushes, warnings shown and warning-filter
    changes, in order (replay() repeats them).

    tee=False (a child): nothing is printed, and the "once" registries start empty, so every
    warning a sequential run could show here is recorded; replay() issues it again where the
    registries are right. tee=True (in-process): everything also happens as usual; the record
    is the same as a child's only if the registries were empty to begin with (fresh())."""

    def __init__(self, tee):
        self.tee, self.events, self.in_warning, self.exact = tee, [], False, True

    def start(self):
        self.out, self.err = sys.stdout, sys.stderr
        self.modules = set(sys.modules)
        self.filters = list(warnings.filters)
        self.show, self.mutated = warnings._showwarnmsg, warnings._filters_mutated
        self.impl, self.showwarning = warnings._showwarnmsg_impl, warnings.showwarning
        if not self.tee:
            for m in list(sys.modules.values()):
                try:
                    reg = vars(m).get("__warningregistry__")
                except TypeError:
                    continue
                if isinstance(reg, dict):
                    reg.clear()
            warnings.onceregistry.clear()
        kind = _Tee if self.tee else _Stream
        sys.stdout, sys.stderr = kind(self, "o", self.out), kind(self, "e", self.err)
        self.hooks = (self._show, self._mutated)  # bound once: each access makes a new object
        warnings._showwarnmsg, warnings._filters_mutated = self.hooks

    def stop(self):
        if (warnings._showwarnmsg, warnings._filters_mutated) != self.hooks \
                or warnings._showwarnmsg_impl is not self.impl or warnings.showwarning is not self.showwarning \
                or list(warnings.filters) != self.filters:
            self.exact = False  # the steps after this one would see other warning settings
        warnings._showwarnmsg, warnings._filters_mutated = self.show, self.mutated
        sys.stdout, sys.stderr = self.out, self.err

    def _mutated(self):
        self.events.append(("m",))
        self.mutated()

    def _show(self, msg):
        if warnings._showwarnmsg_impl is not self.impl or warnings.showwarning is not self.showwarning:
            if warnings.showwarning is not self.showwarning:
                self.exact = False  # a custom display
            return self.show(msg)  # e.g. catch_warnings(record=True) collects it
        module = None
        f = sys._getframe(1)
        while f is not None:  # the frame the warning is attributed to names its module
            if f.f_code.co_filename == msg.filename and f.f_lineno == msg.lineno:
                name = f.f_globals.get("__name__")
                mod = sys.modules.get(name) if isinstance(name, str) else None
                if mod is not None and getattr(mod, "__dict__", None) is f.f_globals and name in self.modules:
                    module = name
                break
            f = f.f_back
        if module is None or msg.file is not None:
            self.exact = False
        self.events.append(("w", msg.message, msg.category, msg.filename, msg.lineno, module))
        if self.tee:
            self.in_warning = True
            try:
                self.show(msg)
            finally:
                self.in_warning = False


def fresh():
    """True when no "once" registry holds a warning the current filters would show, so that a
    step run in-process now shows what it would show in a child (see _Recorder): the registries
    also remember ignored warnings, which change nothing."""
    for m in list(sys.modules.values()):
        try:
            reg = vars(m).get("__warningregistry__")
        except TypeError:
            continue
        if isinstance(reg, dict):
            for key in list(reg):
                if key == "version":
                    continue
                try:
                    text, category, lineno = key
                    if _action(text, category, m.__name__, lineno) != "ignore":
                        return False
                except Exception:
                    return False
    return not warnings.onceregistry


def _action(text, category, module, lineno):
    # the action warnings.warn_explicit takes (a str in a filter must match exactly)
    def matched(pattern, s):
        return pattern is None or (pattern == s if type(pattern) is str else bool(pattern.match(s)))
    for action, msg, cat, mod, ln in warnings.filters:
        if matched(msg, text) and issubclass(category, cat) and matched(mod, module) and (ln == 0 or lineno == ln):
            return action
    return warnings.defaultaction


def _deprecations_hidden():
    """True if every deprecation warning (DeprecationWarning, PendingDeprecationWarning and the
    subclasses defined so far) from library code is ignored, as Python's default filters do.
    Libraries initialize caches once per process, and some warn when they do (matplotlib's
    pyparsing grammars): shown, those warnings would repeat once per child, where a sequential
    run shows them once."""
    classes, todo = set(), [DeprecationWarning, PendingDeprecationWarning]
    while todo:
        c = todo.pop()
        if c not in classes:
            classes.add(c)
            todo.extend(c.__subclasses__())
    for action, msg, cat, mod, ln in warnings.filters:
        main_only = mod == "__main__" or getattr(mod, "pattern", None) in ("__main__", "__main__\\Z")
        if action != "ignore" and not main_only and any(issubclass(c, cat) for c in classes):
            return False
    return all(_action("", c, "a.library", 0) == "ignore" for c in classes)


def record(fn):
    """Runs fn() in-process as usual; returns (result, events), events being what splice()
    needs to repeat its output in a step's child, or None where it cannot be recorded exactly."""
    if not usable() or not fresh():
        return fn(), None
    rec = _Recorder(tee=True)
    rec.start()
    try:
        result = fn()
    finally:
        rec.stop()
    return result, (rec.events if rec.exact else None)


def replay(events):
    """Repeats recorded output: writes, flushes, warnings (through the current filters and
    registries) and warning-filter changes, in order."""
    for ev in events:
        tag = ev[0]
        if tag == "o":
            sys.stdout.write(ev[1])
        elif tag == "e":
            sys.stderr.write(ev[1])
        elif tag == "of":
            sys.stdout.flush()
        elif tag == "ef":
            sys.stderr.flush()
        elif tag == "m":
            warnings._filters_mutated()
        elif tag == "w":
            message, category, filename, lineno, module = ev[1:]
            g = sys.modules[module].__dict__
            warnings.warn_explicit(message, category, filename, lineno, module=module,
                                   registry=g.setdefault("__warningregistry__", {}))
        elif tag == "x":
            replay(ev[1])


def in_step_child():
    return _S["step_child"]


def splice(events):
    """In a step's child: the output recorded by record() is produced here."""
    if not isinstance(sys.stdout, _Stream) or isinstance(sys.stdout, _Tee):
        raise RuntimeError("splice() outside a step's child")
    sys.stdout._rec.events.append(("x", events))


# ---- children --------------------------------------------------------------------------------

def _fork(slot, d=None, parts=None):
    """Forks a child with a raw-output file and a record file: returns (pid, fds) in the parent,
    (0, fds) in the child, or None if creating them failed (the caller then works in-process)."""
    fds = []
    try:
        fds.append(os.memfd_create("sagaconf_raw"))
        fds.append(os.memfd_create("sagaconf_rec"))
        if d is not None:
            os.mkdir(d)
            _S["dirs"].add(d)
        pid = os.fork()
    except (OSError, MemoryError):
        for fd in fds:
            os.close(fd)
        if d is not None and d in _S["dirs"]:
            _S["dirs"].discard(d)
            shutil.rmtree(d, ignore_errors=True)
        _pool().assign(slot, 0)
        return None
    if pid == 0:
        leader = _S["pool"].root == os.getppid()
        _S.update(live={}, dirs=set(), parts=set(), savers=[], step_child=True, depth=0, sig=None, rec=None)
        if leader:
            try:
                os.setpgid(0, 0)  # so that the root can kill the child's whole tree
            except OSError:
                pass
        return 0, fds
    _pool().assign(slot, pid)
    _S["live"][pid] = {"slot": slot, "fds": fds, "status": None, "parts": parts or []}
    if parts:
        _S["parts"].update(t for t, _ in parts)
    return pid, fds


def _child_main(fds, prelude, body):
    """In the child: runs prelude() and body() with output recorded, writes the record, exits."""
    code = 1
    try:
        os.dup2(fds[0], 1)
        os.dup2(fds[0], 2)
        rec = _Recorder(tee=False)
        rec.start()
        _S["rec"] = rec
        try:
            prelude()
            if any(ev[0] != "m" for ev in rec.events):
                rec.exact = False
            del rec.events[:]  # the prelude restores state; its output is the earlier steps'
            ok = body() is not False
        finally:
            rec.stop()
            _S["rec"] = None
        ok = _drain() and ok
        if ok and rec.exact and os.fstat(fds[0]).st_size == 0:
            _write_all(fds[1], pickle.dumps(rec.events, protocol=pickle.HIGHEST_PROTOCOL))
            code = 0
    except BaseException:
        code = 1
    finally:
        try:
            if code:
                _kill_mine()
        finally:
            os._exit(code)


def _write_all(fd, data):
    view = memoryview(data)
    while view:
        view = view[os.write(fd, view):]


def _read_all(fd):
    os.lseek(fd, 0, os.SEEK_SET)
    chunks = []
    while True:
        b = os.read(fd, 1 << 22)
        if not b:
            return b"".join(chunks)
        chunks.append(b)


def _collected(pid, status):
    c = _S["live"][pid]
    c["status"] = status
    _pool().assign(c["slot"], 0)
    c["slot"] = None


def _wait_any():
    """Collects whichever children of this process have exited, sleeping briefly if none has."""
    for _ in range(2):
        got = False
        for pid, c in list(_S["live"].items()):
            if c["status"] is None:
                try:
                    p, status = os.waitpid(pid, os.WNOHANG)
                except ChildProcessError:  # collected elsewhere: its result is unknown
                    p, status = pid, 1 << 8
                if p:
                    _collected(pid, status)
                    got = True
        if got:
            return
        time.sleep(0.01)


def _wait_pid(pid):
    if _S["live"][pid]["status"] is None:
        _collected(pid, os.waitpid(pid, 0)[1])


def _finish(pid):
    """(ok, events) of a collected child; closes its files and forgets it."""
    c = _S["live"].pop(pid)
    ok = os.WIFEXITED(c["status"]) and os.WEXITSTATUS(c["status"]) == 0
    events = None
    if ok:
        try:
            events = pickle.loads(_read_all(c["fds"][1]))
        except Exception:
            ok = False
    for fd in c["fds"]:
        os.close(fd)
    return ok, events


def _descendants(pid):
    found, todo = [], [pid]
    while todo:
        p = todo.pop()
        for line in (_read("/proc/%d/task/%d/children" % (p, p)) or "").split():
            found.append(int(line))
            todo.append(int(line))
    return found


def _kill(pids):
    """Kills these children of this process and everything they forked, collects them and
    removes their partial files."""
    for pid in pids:
        if _S["live"][pid]["status"] is None:
            targets = _descendants(pid)
            if _S["pool"].root == os.getpid():
                try:
                    os.killpg(pid, signal.SIGKILL)  # the child leads its own process group
                except OSError:
                    pass
            for p in [pid] + targets:
                try:
                    os.kill(p, signal.SIGKILL)
                except OSError:
                    pass
    for pid in pids:
        try:
            _wait_pid(pid)
        except ChildProcessError:
            _collected(pid, 1 << 8)
        for t, _ in _S["live"][pid]["parts"]:
            _S["parts"].discard(t)
            try:
                os.remove(t)
            except OSError:
                pass
        _finish(pid)


def _kill_mine():
    """Kills all children of this process and removes their directories and partial files."""
    _kill([pid for pid in _S["live"]])
    _S["savers"] = []
    for d in list(_S["dirs"]):
        shutil.rmtree(d, ignore_errors=True)
    _S["dirs"].clear()
    for p in list(_S["parts"]):
        try:
            os.remove(p)
        except OSError:
            pass
    _S["parts"].clear()


def _on_sigterm(signum, frame):
    prev = _S["sig"]
    _kill_mine()
    signal.signal(signal.SIGTERM, prev if prev is not None else signal.SIG_DFL)
    os.kill(os.getpid(), signal.SIGTERM)


def _enter():
    if _S["depth"] == 0:
        _S["sig"] = None
        try:
            prev = signal.getsignal(signal.SIGTERM)
            if prev is not signal.SIG_IGN and prev is not _on_sigterm:
                signal.signal(signal.SIGTERM, _on_sigterm)
                _S["sig"] = prev
        except (ValueError, OSError):  # not the main thread
            pass
    _S["depth"] += 1


def _leave():
    _S["depth"] -= 1
    if _S["depth"] == 0 and _S["sig"] is not None:
        try:
            signal.signal(signal.SIGTERM, _S["sig"])
        except (ValueError, OSError):
            pass
        _S["sig"] = None


def _merge(src, dst):
    for name in sorted(os.listdir(src)):
        s, t = os.path.join(src, name), os.path.join(dst, name)
        if os.path.isdir(t) and not os.path.islink(t) and os.path.isdir(s) and not os.path.islink(s):
            _merge(s, t)
        else:
            os.replace(s, t)
    os.rmdir(src)


def run_steps(steps, savedir, preludes, costs, mem, inprocess=(), restart=None):
    """Same as `for step in steps: step(savedir)`, run in parallel.

    preludes[i]() restores, in a process that has run none of the steps, the global state that
    steps[:i] leave behind (preludes has one more entry, for after the last step). costs[i]
    orders the launches (largest first); mem[i] is step i's estimated peak extra memory in a
    child. Steps listed in `inprocess` always run in this process, in their turn. A step that
    returns False (a failure it reported) or fails in its child is the last step whose children
    are kept: the steps after it run in-process. restart[j] = k (k < j): if step j fails in its
    child, the run is the sequential one from step k on, not from step j.
    """
    if not usable() or len(steps) < 2:
        for step in steps:
            step(savedir)
        return
    _enter()
    mine = {"done": {}}
    try:
        _run(steps, savedir, preludes, costs, mem, set(inprocess), dict(restart or {}), mine)
    except BaseException:
        if "tag" in mine:
            _discard(mine, 0, savedir, mine["done"])
        raise
    finally:
        _leave()


def _step_dir(savedir, tag, i):
    return os.path.join(savedir, ".sagaconf_step%s_%d" % (tag, i))


def _run(steps, savedir, preludes, costs, mem, inprocess, restart, mine):
    n = len(steps)
    pool = _pool()
    pending = sorted((i for i in range(n) if i not in inprocess), key=lambda i: -costs[i])
    nofork, tried = False, 0.0
    _S["seq"] += 1
    mine["tag"] = tag = "%d_%d" % (os.getpid(), _S["seq"])
    pids = mine.setdefault("pids", {})
    done = mine["done"]
    i = 0
    while i < n:
        if pending and not nofork and time.monotonic() - tried > 0.2:
            tried = time.monotonic()
            for j in list(pending):  # launch what fits, most costly first
                slot = pool.reserve(mem[j])
                if slot is None:
                    continue
                d = _step_dir(savedir, tag, j)
                got = _fork(slot, d)
                if got is None:
                    nofork = True
                    break
                pending.remove(j)
                pid, fds = got
                if pid == 0:
                    _child_main(fds, preludes[j], lambda: steps[j](d))
                pids[j] = pid
                tried = 0.0

        if i in pids or i in done:  # commit step i, once it and the steps restarting from it are done
            waiting = False
            for j in [i] + [j for j in range(i + 1, n) if restart.get(j, j) <= i]:
                if j in done:
                    continue
                if j in pids:
                    if _S["live"][pids[j]]["status"] is None:
                        waiting = True
                    else:
                        done[j] = _finish(pids.pop(j))
                elif j in pending and not nofork and _running():
                    waiting = True
                else:  # step j will run in-process, so step i must too
                    done[j] = (False, None)
            if waiting:
                _wait_any()
                continue
            if any(not done[j][0] for j in done if j >= i):
                _discard(mine, i, savedir, done)
                preludes[i]()
                for step in steps[i:]:
                    step(savedir)
                return
            replay(done.pop(i)[1])
            d = _step_dir(savedir, tag, i)
            _merge(d, savedir)
            _S["dirs"].discard(d)
            i += 1
            continue

        # step i was not forked: it runs here, unless a running child may free what it needs
        if i not in inprocess and not nofork and _running():
            _wait_any()
            continue
        if i in pending:
            pending.remove(i)
        preludes[i]()
        i += 1
        if steps[i - 1](savedir) is False:
            _discard(mine, i, savedir, done)
            for step in steps[i:]:
                step(savedir)
            return
    preludes[n]()


def _running():
    return any(c["status"] is None for c in _S["live"].values())


def _discard(mine, first, savedir, done=None):
    """Kills the children of steps first.. of one run_steps call (forgetting those collected
    already) and removes their directories."""
    pids = mine.get("pids", {})
    _kill([pids[j] for j in sorted(pids) if j >= first and pids[j] in _S["live"]])
    doomed = [j for j in pids if j >= first] + [j for j in (done or {}) if j >= first]
    for j in doomed:
        pids.pop(j, None)
        if done:
            done.pop(j, None)
        d = _step_dir(savedir, mine["tag"], j)
        _S["dirs"].discard(d)
        shutil.rmtree(d, ignore_errors=True)


@contextlib.contextmanager
def scratch_file(path):
    """Yields path, a file to work in that is removed afterwards, also when this process is
    interrupted or stopped with SIGTERM while running steps, or when its step is discarded."""
    _S["parts"].add(path)
    try:
        yield path
    finally:
        _S["parts"].discard(path)
        try:
            os.remove(path)
        except FileNotFoundError:
            pass


# ---- figures ---------------------------------------------------------------------------------

def savefig_pdf_svg(pdf, svg):
    """Same as plt.savefig(pdf, format='pdf'); plt.savefig(svg, format='svg').

    In a step's child, a forked child does both while the step goes on; its record goes where
    the step called this, and the step's child waits for it before it ends. Elsewhere, a forked
    child saves the SVG from the same figure state while this process saves the PDF; it must
    print and warn nothing, else the SVG is saved again here. (Saving draws the figure without
    changing what a later save or plot draws.) If a forked save fails, the SVG is saved here,
    or, in a step's child, the step fails and runs again in-process.
    """
    from matplotlib import pyplot as plt
    parts = [(pdf + ".sagaconf_part", pdf), (svg + ".sagaconf_part", svg)]
    if usable() and not any(os.path.lexists(t) for t, _ in parts):
        if _S["step_child"]:
            if _S["rec"] is not None and _save_async(plt, parts):
                return
        elif _save_svg_forked(plt, pdf, parts[1]):
            return
    plt.savefig(pdf, format='pdf')
    plt.savefig(svg, format='svg')


def _save_async(plt, parts):
    while sum(_S["live"][p]["status"] is None for p in _S["savers"]) >= MAX_SAVERS:
        _wait_any()
    slot = _pool().reserve(SAVER_BYTES)
    if slot is None:
        return False
    got = _fork(slot, parts=parts)
    if got is None:
        return False
    pid, fds = got
    if pid == 0:
        def body():
            plt.savefig(parts[0][0], format='pdf')
            plt.savefig(parts[1][0], format='svg')
        _child_main(fds, lambda: None, body)
    place = []  # the save's record, filled in by _drain
    _S["rec"].events.append(("x", place))
    _S["live"][pid]["place"] = place
    _S["savers"].append(pid)
    return True


def _drain():
    """Waits for this step child's figure saves; if they all succeeded, puts their records in
    place and renames their files into place (else removes the files). True if they did."""
    ok = True
    for pid in _S["savers"]:
        _wait_pid(pid)
    for pid in _S["savers"]:
        parts, place = _S["live"][pid]["parts"], _S["live"][pid]["place"]
        good, events = _finish(pid)
        ok = ok and good
        if good:
            place.extend(events)
        for t, final in parts:
            _S["parts"].discard(t)
            if ok:
                os.replace(t, final)
            elif os.path.lexists(t):
                os.remove(t)
    _S["savers"] = []
    return ok


def _save_svg_forked(plt, pdf, part):
    slot = _pool().reserve(SAVER_BYTES)
    if slot is None:
        return False
    _enter()
    try:
        got = _fork(slot, parts=[part])
        if got is None:
            return False
        pid, fds = got
        if pid == 0:
            _child_main(fds, lambda: None, lambda: plt.savefig(part[0], format='svg'))
        try:
            plt.savefig(pdf, format='pdf')
            _wait_pid(pid)
        except BaseException:
            _kill([pid])
            raise
    finally:
        _leave()
    ok, events = _finish(pid)
    _S["parts"].discard(part[0])
    if ok and not events:
        os.replace(part[0], part[1])
    else:
        if os.path.lexists(part[0]):
            os.remove(part[0])
        plt.savefig(part[1], format='svg')
    return True
