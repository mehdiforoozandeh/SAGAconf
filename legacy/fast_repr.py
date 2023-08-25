"""
call from terminal to get reproducibility results for any pair of annotations

gets the following arguments:

TYPE -> short OR long
ANNOTATION 1 DIR -> <dir>
ANNOTATION 2 DIR -> <dir>
OUTDIR -> <dir>

example:

    python fast_repr.py short segway_runs/GM12878_rep1 segway_runs/GM12878_rep2 repr_GM12878_rep1vsrep2
"""

from run import *
import sys

input_dict = {
    "runtype": sys.argv[1],
    "rep1_dir": sys.argv[2],
    "rep2_dir": sys.argv[3],
    "output_dir": sys.argv[4]
}

if os.path.exists(input_dict["output_dir"]) == False:
    os.mkdir(input_dict["output_dir"])

if input_dict["runtype"] == "short":
    with open(input_dict["output_dir"]+"/run_info.txt", 'w') as fw:
        fw.write(str(input_dict))

    if "segway" in input_dict["rep1_dir"]:
        type = "segw"
    elif "chromhmm" in input_dict["rep1_dir"]:
        type = "chmm"

    get_short_report(input_dict["rep1_dir"], input_dict["rep2_dir"], input_dict["output_dir"], type=type, mnemons=False)


elif input_dict["runtype"] == "long":
    with open(input_dict["output_dir"]+"/run_info.txt", 'w') as fw:
        fw.write(str(input_dict))
    full_reproducibility_report(input_dict["rep1_dir"], input_dict["rep2_dir"], input_dict["output_dir"], run_on_subset=False, mnemons=False)

else:
    print("did not specify correct run type (short or long)")