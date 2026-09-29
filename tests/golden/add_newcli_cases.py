"""Add a `<name>_newcli` twin to every case: the same run through the new `sagaconf` CLI
syntax, which must match the legacy script byte for byte.

    python tests/golden/add_newcli_cases.py tests/golden/cases.json
"""
import json
import sys

MODE_FLAGS = {"-q": "quick", "--quick": "quick", "--r_only": "rvalues", "--ct_only": "celltype",
              "--v_seglength": "seglength", "--merge_only": "merge", "--active_regions": "active-regions"}
RENAME = {"-bm": "--base-mnemonics", "-vm": "--verif-mnemonics", "-s": "--chr21-only",
          "-w": "--window-size", "-to": "--iou-threshold", "-tr": "--repr-threshold",
          "-k": "--merge-k", "-v": "--verbose", "--out_format": "--out-format"}


def translate(case):
    args, out, modes = case["args"], [], []
    for a in args:
        if a in MODE_FLAGS:
            modes.append(MODE_FLAGS[a])
        else:
            out.append(RENAME.get(a, a))
    if case["script"] == "SAGAconf_parser.py":
        return ["parse"] + out
    if "active-regions" in modes:
        # The legacy script reads these cwd-relative paths; pass the same ones explicitly.
        out = ["--ccre-file", "src/biointerpret/GRCh38-cCREs.bed",
               "--meuleman-file", "src/biointerpret/Meuleman.tsv"] + out
    return ["run"] + (["--mode", modes[0]] if modes else []) + out


def main(path):
    cases = [c for c in json.load(open(path)) if not c["name"].endswith("_newcli")]
    out = []
    for c in cases:
        out.append(c)
        twin = dict(c, name=c["name"] + "_newcli", new_script="sagaconf/__main__.py", new_args=translate(c))
        out.append(twin)
    json.dump(out, open(path, "w"), indent=1)
    print("%s: %d cases" % (path, len(out)))


if __name__ == "__main__":
    main(sys.argv[1])
