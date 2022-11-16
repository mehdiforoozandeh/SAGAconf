"""
source  ~/miniconda3/etc/profile.d/conda.sh
conda activate p27segenv
python _interpret.py sigdist
conda deactivate
"""
import os

def CROSS_chmm_sigdist(chmmruns_dir, original_files_dir):
    ls0 = os.listdir(chmmruns_dir)
    ls1 = os.listdir(original_files_dir)

    for ct in ls1:
        for run in ls0:
            try:
                if ct in run:
                    ls2 = os.listdir("{}/{}/".format(chmmruns_dir, run))

                    if "concat" in run:
                        for ff in ls2:
                            if "dense.bed" in ff:
                                if "rep1" in ff and "psd" not in ff:
                                    segbed = "{}/{}/{}".format(chmmruns_dir, run, ff)
                                    gd = "{}/{}/rep2.genomedata".format(original_files_dir, ct)
                                    outdir = "{}/{}/{}".format(chmmruns_dir, run, 'CROSS_sigdist_rep1')
                                    if os.path.exists(outdir) == False:
                                        os.system(
                                            'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, outdir))

                                elif "rep2" in ff and "psd" not in ff:
                                    segbed = "{}/{}/{}".format(chmmruns_dir, run, ff)
                                    gd = "{}/{}/rep1.genomedata".format(original_files_dir, ct)
                                    outdir = "{}/{}/{}".format(chmmruns_dir, run, 'CROSS_sigdist_rep2')
                                    if os.path.exists(outdir) == False:
                                        os.system(
                                            'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, outdir))

                    else:
                        for ff in ls2:
                            if "dense.bed" in ff:
                                if "rep1" in ff and "psd" not in ff:
                                    segbed = "{}/{}/{}".format(chmmruns_dir, run, ff)
                                    gd = "{}/{}/rep2.genomedata".format(original_files_dir, ct)
                                    outdir = "{}/{}/{}".format(chmmruns_dir, run, 'CROSS_sigdist')
                                    if os.path.exists(outdir) == False:
                                        os.system(
                                            'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, outdir))

                                elif "rep2" in ff and "psd" not in ff:
                                    segbed = "{}/{}/{}".format(chmmruns_dir, run, ff)
                                    gd = "{}/{}/rep1.genomedata".format(original_files_dir, ct)
                                    outdir = "{}/{}/{}".format(chmmruns_dir, run, 'CROSS_sigdist')
                                    if os.path.exists(outdir) == False:  
                                        os.system(
                                            'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, outdir))
            except:
                pass

def CROSS_segway_sigdist(segwayruns_dir, originalfiles_dir): 
    ls0 = os.listdir(originalfiles_dir)
    ls1 = os.listdir(segwayruns_dir)
    for i in ls0:
        for j in ls1:
            # try:
            segbed = segwayruns_dir+j+'/segway.bed'
            if i in j:
                if "concat" in j:

                    if "rep1" in j and "psd" not in j:
                        gd = originalfiles_dir+i+"/concat_rep2.genomedata"
                    elif "rep2" in j and "psd" not in j:
                        gd = originalfiles_dir+i+"/concat_rep1.genomedata"
                    
                    print(gd)
                    
                    if os.path.exists(segwayruns_dir+j+'/sigdist') == False:
                        os.system(
                            'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, segwayruns_dir+j+'/CROSS_sigdist'))

                else:
                    if "rep1" in j and "psd" not in j:
                        gd = originalfiles_dir+i+"/rep2.genomedata"
                    elif "rep2" in j and "psd" not in j:
                        gd = originalfiles_dir+i+"/rep1.genomedata"

                    print(gd)
                    
                    if os.path.exists(segwayruns_dir+j+'/sigdist') == False:
                        os.system(
                            'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, segwayruns_dir+j+'/CROSS_sigdist'))
            # except:
            #     pass


if __name__=="__main__":
    originalfiles = "files/"
    chmmruns = "chromhmm_runs/"
    segwayruns = "segway_runs/"

    # print("running chromhmm signal dist")
    # CROSS_chmm_sigdist(chmmruns_dir=chmmruns, original_files_dir=originalfiles)
        
    print("running segway signal dist")
    CROSS_segway_sigdist(segwayruns, originalfiles)
