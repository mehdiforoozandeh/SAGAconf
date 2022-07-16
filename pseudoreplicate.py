import os, shutil, sys, re
import pandas as pd
from sklearn.utils import shuffle
import requests, os, itertools, ast


"""
STARTING FROM BAM FILES

- bam to bed
- two random subsamples from the initial bed
    - randomly [shuffle] and split reads
- save new pseudoreplicates on a new bed
- bed to fcoc
- use bed for chmm
- use fcoc for segway

"""

def bam_to_bed(bamfile):
    os.system("bedtools bamtobed -i {} > {}".format(
        bamfile, bamfile.replace(".bam", ".bed")
    ))
    
def read_bed(bedfile):
    bed = pd.read_csv(bedfile, sep="\t", header=None)
    bed.columns=["chr", "start", "end", "-", "MAPQ", "strand"]
    return bed

def sample_pseudo(bed, ignore_decoy=True):
    if ignore_decoy:
        bed = bed.drop(bed[bed.chr == "chrEBV"].index).reset_index(drop=True)

    bed = shuffle(bed, random_state=5).reset_index(drop=True)
    splitting_point = int(len(bed)/2)
    bed1 = bed.iloc[:splitting_point,:].sort_values(by=["chr","start"]).reset_index(drop=True)
    bed2 = bed.iloc[splitting_point:,:].sort_values(by=["chr","start"]).reset_index(drop=True)
    return bed1, bed2

def save_pseudo(psd_bed1, psd_bed2, initial_bam_dir):
    psd_bed1.to_csv(initial_bam_dir.replace(".bam", "_psdrep1.bed"), index=False, header=False, sep="\t")
    psd_bed2.to_csv(initial_bam_dir.replace(".bam", "_psdrep2.bed"), index=False, header=False, sep="\t")

def genomesize(chrsz):
    chr = pd.read_csv(chrsz, sep="\t", header=None)
    return chr.iloc[:,1].sum()

def bed_to_fc(bedfile, fraglen, chrsz, gensz, outdir):
    os.system(
        "python bam_to_fc.py {} --fraglen {} --shift 0 --chrsz {} --gensz {} --pval-thresh 0.01 --ctl-subsample 1 --out-dir {} --log-level INFO".format(
            bedfile, fraglen, chrsz, gensz, outdir
        )
    )

def bw2bg(dir):
    ls_bws = [dir+"/"+e for e in os.listdir(dir) if ".bigwig" in e]
    for e in ls_bws:
        os.system(
            "bigWigToBedGraph {} {}".format(
                e, e.replace(".bed.",".").replace(".bigwig", ".bedGraph")
            )
        )

def clean_up(dir):
    ls_bgs = [dir+"/"+e for e in os.listdir(dir) if ".bedGraph" in e]
    parent_dir = "/".join(dir.split("/")[:-1]) + "/"
    for e in ls_bgs:
        os.system(
            "mv {} {}".format(
                e, parent_dir
            )
        )
    
    shutil.rmtree(dir)

def get_fraglen(encode_file_accession):
    """
    searches encode json outputs to find out the fragment length or
    estimated fragment length corresponding to each bam file.
    """
    headers = {'accept': 'application/json'}

    base_search_url = """https://www.encodeproject.org/files/{}/?format=json""".format(encode_file_accession)
    
    # GET the object
    response = requests.get(base_search_url, headers=headers)

    # Extract the JSON response as a Python dictionary
    search_results = response.json()
    fraglen = str(search_results["quality_metrics"])[
            str(search_results["quality_metrics"]).find("fragment_len"): 
            str(search_results["quality_metrics"]).find("fragment_len")+18].split(":")[1]

    print(
        encode_file_accession, int(fraglen))
    return int(fraglen)
    

def psdrep_pipeline(initial_bam):
    accession = initial_bam.split("/")[-1].replace(".bam", "")

    if os.path.exists(initial_bam.replace(".bam", ".bed")) == False:
        print("converting bam to bed ...")
        bam_to_bed(initial_bam)

    print("subsampling pseudoreplicates ...")
    psd_bed1, psd_bed2 = sample_pseudo(read_bed(initial_bam.replace(".bam", ".bed")))

    fraglen = int(get_fraglen(accession))

    gensz = genomesize("hg38.txt")

    print("fragment length: {}".format(fraglen))
    print("genome size: {}".format(gensz))

    print("saving psedoreplicates...")
    save_pseudo(psd_bed1, psd_bed2, initial_bam)
    del psd_bed1, psd_bed2

    print("getting signals for psdrep1")
    bed_to_fc(
        initial_bam.replace(".bam", "_psdrep1.bed"), 
        fraglen, "hg38.txt", gensz,
        initial_bam.replace(".bam", "psdrep1_signals"))

    print("getting signals for psdrep2")
    bed_to_fc(
        initial_bam.replace(".bam", "_psdrep2.bed"), 
        fraglen, "hg38.txt", gensz,
        initial_bam.replace(".bam", "psdrep2_signals"))

    print("converting bigwig to bedgraph and Cleaning up")
    bw2bg(initial_bam.replace(".bam", "psdrep1_signals"))
    clean_up(initial_bam.replace(".bam", "psdrep1_signals"))

    bw2bg(initial_bam.replace(".bam", "psdrep2_signals"))
    clean_up(initial_bam.replace(".bam", "psdrep2_signals"))


# if __name__=="__main__":
#     acclist = ['ENCFF891JRY', 'ENCFF760VTC', 'ENCFF992JKA', 'ENCFF766LEF', 'ENCFF014BCD', 'ENCFF616FAT', 'ENCFF872QFO', 'ENCFF656YCM', 'ENCFF630ANH', 'ENCFF113FCH', 'ENCFF785GAD', 'ENCFF596XSW', 'ENCFF150SGV', 'ENCFF183EDU', 'ENCFF629ODR', 'ENCFF520XFF', 'ENCFF141KHX', 'ENCFF278FTJ', 'ENCFF573MNS', 'ENCFF310NWI', 'ENCFF384IVB', 'ENCFF007YZT', 'ENCFF591CTO', 'ENCFF874SMO', 'ENCFF144PKE', 'ENCFF121RHF', 'ENCFF640GCK', 'ENCFF907MNY', 'ENCFF159LCH', 'ENCFF392ZKG', 'ENCFF463WUU', 'ENCFF905CZD', 'ENCFF326KCT', 'ENCFF272JVI', 'ENCFF775XVG', 'ENCFF880HKV', 'ENCFF168UCU', 'ENCFF352HXD', 'ENCFF537BKE', 'ENCFF415GHS', 'ENCFF238ZAC', 'ENCFF446FUS', 'ENCFF899MAF', 'ENCFF583GSH', 'ENCFF256UUQ', 'ENCFF564SVK', 'ENCFF703AIM', 'ENCFF685PPQ', 'ENCFF201OHW', 'ENCFF465UWC', 'ENCFF269GKF', 'ENCFF863HNL', 'ENCFF565DCK', 'ENCFF149MXA', 'ENCFF633BHN', 'ENCFF505NMT', 'ENCFF353YPB', 'ENCFF104THG', 'ENCFF677MAG', 'ENCFF155UQU', 'ENCFF047NLO', 'ENCFF255QRL', 'ENCFF385FLM', 'ENCFF923YUN', 'ENCFF008LEY', 'ENCFF597OEX', 'ENCFF081ODV', 'ENCFF687JRM', 'ENCFF126HII', 'ENCFF183HME', 'ENCFF843BWY', 'ENCFF626ORM', 'ENCFF465NXQ', 'ENCFF391BOB', 'ENCFF967CHL', 'ENCFF785NUK', 'ENCFF405REQ', 'ENCFF413QYQ', 'ENCFF729HQN', 'ENCFF744JQU', 'ENCFF306ENK', 'ENCFF551PNK', 'ENCFF889UPU', 'ENCFF747AWB', 'ENCFF132EDT', 'ENCFF592EVS', 'ENCFF334NVA', 'ENCFF748ISL', 'ENCFF650RBB', 'ENCFF134YPW', 'ENCFF822OCY', 'ENCFF680ERD', 'ENCFF609ZAE', 'ENCFF101TRI', 'ENCFF711QAI', 'ENCFF374GQR', 'ENCFF299WXR', 'ENCFF312WJY', 'ENCFF620CCL', 'ENCFF654NDM', 'ENCFF008RRT', 'ENCFF403SOH', 'ENCFF317FFS', 'ENCFF571MVG', 'ENCFF712AAP', 'ENCFF249TCH', 'ENCFF826OLG', 'ENCFF418SFP', 'ENCFF622IIB', 'ENCFF909MGB', 'ENCFF706FSX', 'ENCFF923NDP', 'ENCFF300IWJ', 'ENCFF491EQT', 'ENCFF650IXI', 'ENCFF741HOA']
#     for acc in acclist:
#         try:
#             get_fraglen(acc)
#         except:
#             print(acc)
