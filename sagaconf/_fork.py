"""Run independent report steps in forked children, with the output of a sequential run.

Each step writes into its own directory under savedir and prints into its own in-memory
files. Once every child has finished, the steps are committed in step order: printed text is
replayed to stdout/stderr and files are moved into savedir. A step that fails in its child is
run again here, in-process, after the global state of the steps before it is restored, and so
are the steps after it: from the first failure on, the run is the sequential one.
"""
import os
import shutil
import sys
import time

MAX_WORKERS = 4


def _workers():
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


def _usable():
    if not hasattr(os, "fork") or not hasattr(os, "memfd_create") or _workers() < 2:
        return False
    try:
        return (sys.stdout is sys.__stdout__ and sys.stderr is sys.__stderr__
                and sys.stdout.fileno() == 1 and sys.stderr.fileno() == 2)
    except Exception:
        return False


def _child(prelude, step, d, fds):
    code = 0
    try:
        os.dup2(fds[0], 1)
        os.dup2(fds[1], 2)
        prelude()
        _async["pid"], _async["pending"] = os.getpid(), []
        if step(d) is False:  # a step that handles its own failure reports it this way
            code = 1
    except BaseException:
        code = 1
    try:
        sys.stdout.flush()
        sys.stderr.flush()
    except BaseException:
        code = 1
    if not _drain():
        code = 1
    os._exit(code)


def _merge(src, dst):
    for name in sorted(os.listdir(src)):
        s, t = os.path.join(src, name), os.path.join(dst, name)
        if os.path.isdir(t) and not os.path.islink(t) and os.path.isdir(s):
            _merge(s, t)
        else:
            os.replace(s, t)
    os.rmdir(src)


def _replay(fd, stream):
    os.lseek(fd, 0, os.SEEK_SET)
    chunks = []
    while True:
        b = os.read(fd, 1 << 20)
        if not b:
            break
        chunks.append(b)
    os.close(fd)
    if chunks:
        stream.flush()
        stream.buffer.write(b"".join(chunks))
        stream.buffer.flush()


def run_steps(steps, savedir, preludes, costs, workers=MAX_WORKERS, quiet=False):
    """Same as `for i, step in enumerate(steps): step(savedir)`, run in parallel.

    preludes[i]() restores, in a process that has run none of the steps, the global state
    steps[:i] leave behind (preludes has one more entry, for after the last step); costs[i]
    orders the launches (largest first). With quiet=True, a step that writes to stderr (a
    warning, which a sequential run may print only once) counts as failed.
    """
    if not _usable() or len(steps) < 2:
        for step in steps:
            step(savedir)
        return

    dirs = [os.path.join(savedir, ".sagaconf_step%d" % i) for i in range(len(steps))]
    if any(os.path.lexists(d) for d in dirs):
        for step in steps:
            step(savedir)
        return

    sys.stdout.flush()
    sys.stderr.flush()
    pending = sorted(range(len(steps)), key=lambda i: -costs[i])
    workers = min(workers, _workers())
    running, fds, ok = {}, {}, {}
    while pending or running:
        while pending and len(running) < workers:
            i = pending.pop(0)
            fds[i] = (os.memfd_create("sagaconf_out"), os.memfd_create("sagaconf_err"))
            os.mkdir(dirs[i])
            pid = os.fork()
            if pid == 0:
                _child(preludes[i], steps[i], dirs[i], fds[i])
            running[pid] = i
        done = False
        for pid in list(running):
            p, status = os.waitpid(pid, os.WNOHANG)
            if p:
                i = running.pop(pid)
                ok[i] = (os.WIFEXITED(status) and os.WEXITSTATUS(status) == 0
                         and not (quiet and os.fstat(fds[i][1]).st_size))
                done = True
        if not done:
            time.sleep(0.02)

    for i in range(len(steps)):
        if not ok[i]:
            for j in range(i, len(steps)):
                os.close(fds[j][0])
                os.close(fds[j][1])
                shutil.rmtree(dirs[j], ignore_errors=True)
            preludes[i]()
            for step in steps[i:]:
                step(savedir)
            return
        _replay(fds[i][0], sys.stdout)
        _replay(fds[i][1], sys.stderr)
        _merge(dirs[i], savedir)
    preludes[len(steps)]()


# In a step's child, figures are saved by their own children while the step goes on; the step's
# child waits for them (and renames their files into place) before it finishes.
_async = {"pid": None, "pending": []}


def _drain():
    ok = True
    for pid, fds, parts in _async["pending"]:
        ok = _svg_child_ok(pid, fds) and ok
    for pid, fds, parts in _async["pending"]:
        for tmp, final in parts:
            if ok:
                os.replace(tmp, final)
            elif os.path.lexists(tmp):
                os.remove(tmp)
    _async["pending"] = []
    return ok


def savefig_pdf_svg(pdf, svg):
    """Same as plt.savefig(pdf, format='pdf'); plt.savefig(svg, format='svg')."""
    parts = [(pdf + ".sagaconf_part", pdf), (svg + ".sagaconf_part", svg)]
    if _async["pid"] != os.getpid() or not _usable() or any(os.path.lexists(t) for t, _ in parts):
        return _savefig_pdf_svg(pdf, svg)
    from matplotlib import pyplot as plt
    sys.stdout.flush()
    sys.stderr.flush()
    fds = (os.memfd_create("sagaconf_out"), os.memfd_create("sagaconf_err"))
    pid = os.fork()
    if pid == 0:
        code = 0
        try:
            os.dup2(fds[0], 1)
            os.dup2(fds[1], 2)
            _async["pid"] = None
            _savefig_pdf_svg(parts[0][0], parts[1][0], last=True)
            sys.stdout.flush()
            sys.stderr.flush()
        except BaseException:
            code = 1
        os._exit(code)
    _async["pending"].append((pid, fds, parts))


def _savefig_pdf_svg(pdf, svg, last=False):
    """Same as plt.savefig(pdf, format='pdf'); plt.savefig(svg, format='svg').

    Saving draws the figure without changing it (no layout engine is left set), so the SVG is
    drawn in a forked child, from the same figure state, into a temporary file while this
    process draws the PDF; it is renamed to svg once the PDF is saved. If the child fails or
    prints anything, the SVG is saved here, after the PDF, instead.
    """
    from matplotlib import pyplot as plt
    pid, tmp = None, svg + ".sagaconf_part"
    if last:
        _skip_predraw(plt.gcf())
    if _usable() and not os.path.lexists(tmp):
        sys.stdout.flush()
        sys.stderr.flush()
        fds = (os.memfd_create("sagaconf_out"), os.memfd_create("sagaconf_err"))
        pid = os.fork()
        if pid == 0:
            code = 0
            try:
                os.dup2(fds[0], 1)
                os.dup2(fds[1], 2)
                _skip_predraw(plt.gcf())
                plt.gcf().savefig(tmp, format='svg')  # this process ends after saving
                sys.stdout.flush()
                sys.stderr.flush()
            except BaseException:
                code = 1
            os._exit(code)
    try:
        if last:  # the process ends after saving: pyplot's redraw after saving changes nothing
            plt.gcf().savefig(pdf, format='pdf')
        else:
            plt.savefig(pdf, format='pdf')
    finally:
        ok = pid is not None and _svg_child_ok(pid, fds)
        if pid is not None and (not ok or sys.exc_info()[0] is not None) and os.path.lexists(tmp):
            os.remove(tmp)
    if ok:
        os.replace(tmp, svg)
    else:
        plt.savefig(svg, format='svg')


def _skip_predraw(fig):
    # tight_layout leaves a PlaceHolderLayoutEngine, which makes savefig draw the figure once
    # with drawing disabled before drawing it for real; that pass lays out nothing (the
    # placeholder's execute is a no-op) and the real draw does the same work again. In a
    # process that ends after saving, the figure is saved without it.
    from matplotlib.layout_engine import PlaceHolderLayoutEngine
    if type(fig.get_layout_engine()) is PlaceHolderLayoutEngine:
        fig.set_layout_engine(None)


def _svg_child_ok(pid, fds):
    status = os.waitpid(pid, 0)[1]
    ok = (os.WIFEXITED(status) and os.WEXITSTATUS(status) == 0
          and not os.fstat(fds[0]).st_size and not os.fstat(fds[1]).st_size)
    os.close(fds[0])
    os.close(fds[1])
    return ok
