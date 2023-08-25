import pstats
import requests, os, sys

def gtf_file(gtf_filename = 'biointerpret/gencode.v29.primary_assembly.annotation_UCSC_names.gtf'): 
    if not os.path.isfile(gtf_filename):
        url='https://www.encodeproject.org/files/\
        gencode.v29.primary_assembly.annotation_UCSC_names/\
            @@download/gencode.v29.primary_assembly.annotation_UCSC_names.gtf.gz'

        r = requests.get(url, allow_redirects=True)
        open(gtf_filename + '.gz', 'wb').write(r.content)

        os.system('gzip -d {}'.format(gtf_filename+ '.gz'))

        return gtf_filename

def __feature_aggreg(exp_name, segbed, gtf):
    '''
    use segtools to run segtools feat_aggreg
    '''
    outdir = 'label_interpretation/feature_files/{}'.format(exp_name)

    os.system('segtools-aggregation --normalize --mode=gene {} {} --outdir={}'.format(segbed, gtf, outdir))

def __signal_dist(exp_name, segbed, gd):
    '''
    use segtools to run segtools signal_distribution
    '''
    outdir = 'label_interpretation/feature_files/{}'.format(exp_name)

    os.system('segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, outdir))

def segway_concat_beds(segwayruns_dir):
    "takes care of the added index infront of the chromosomes in the bed file"
    ls0 = os.listdir(segwayruns_dir)

    for l in ls0:
        if "concat" in l:
            ls1 = os.listdir(segwayruns_dir+"/"+l+"/")
            if "segway_temp.bed" not in ls1 and "segway.bed" in ls1:
                
                with open(segwayruns_dir+"/"+l+"/segway.bed", 'r') as fr:
                    lines = fr.readlines()
                    with open(segwayruns_dir+"/"+l+"/segway_temp.bed", 'w') as fw:
                        for ll in lines:
                            if ll[:3] == "chr":
                                lll = ll.split("\t")
                                lll[0] = lll[0][:-2]
                                lll = "\t".join(lll)
                                fw.write(lll)

                            else:
                                fw.write(ll)

def segway_sigdist(segwayruns_dir, originalfiles_dir): 
    ls0 = os.listdir(originalfiles_dir)
    ls1 = os.listdir(segwayruns_dir)
    for i in ls0:
        for j in ls1:
            segbed = segwayruns_dir+j+'/segway.bed'
            if i in j:
                if "concat" in j:

                    if "rep1" in j:
                        gd = originalfiles_dir+i+"/concat_rep1.genomedata"
                    elif "rep2" in j:
                        gd = originalfiles_dir+i+"/concat_rep2.genomedata"

                elif "_rep1_psdrep1" in j:
                    gd = originalfiles_dir+i+"/rep1_psdrep1.genomedata"
                elif "_rep1_psdrep2" in j:
                    gd = originalfiles_dir+i+"/rep1_psdrep2.genomedata"
                elif "_rep2_psdrep1" in j:
                    gd = originalfiles_dir+i+"/rep2_psdrep1.genomedata"
                elif "_rep2_psdrep2" in j:
                    gd = originalfiles_dir+i+"/rep2_psdrep2.genomedata"

                else:
                    if "rep1" in j:
                        gd = originalfiles_dir+i+"/rep1.genomedata"
                    elif "rep2" in j:
                        gd = originalfiles_dir+i+"/rep2.genomedata"

                if os.path.exists(segwayruns_dir+j+'/sigdist') == False:
                    os.system(
                        'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, segwayruns_dir+j+'/sigdist'))

def segway_feataggr(
    segwayruns_dir, gtffile="biointerpret/gencode.v29.primary_assembly.annotation_UCSC_names.gtf"):
    ls0 = os.listdir(segwayruns_dir)
    for i in ls0:
        if "concat" in i:
            segbed = segwayruns_dir+i+'/segway_temp.bed'
        else:
            segbed = segwayruns_dir+i+'/segway.bed'
            
        if os.path.exists("{}/{}/aggre/feature_aggregation.tab".format(segwayruns_dir, i)) == False:
            os.system(
                'segtools-aggregation --normalize --mode=gene {} {} --outdir={}'.format(
                    segbed, gtffile, "{}/{}/aggre".format(segwayruns_dir, i)))

def segway_get_mnem(segwayruns_dir):
    """
    for each run:
        mkdir(feature_files/runname)
        cp run/aggre/.tab feature_files/runname/.tab
        cp run/sigdist/.tab feature_files/runname/.tab
    
    run apply.py mnems

    for each run:
        cp mnems/run/mnem.txt run/mnem.txt
        """
    ls1 = os.listdir(segwayruns_dir)
    for run in ls1:
        if os.path.exists("{}/{}/sigdist/signal_distribution.tab".format(segwayruns_dir, run)) and os.path.exists("{}/{}/aggre/feature_aggregation.tab".format(segwayruns_dir, run)):
            with open("{}/{}/sigdist/signal_distribution.tab".format(segwayruns_dir, run), 'r') as file:
                lines = file.readlines()
                lines = "".join(lines)
            if "nan" not in lines:
                if os.path.exists("biointerpret/feature_files/{}".format(run))==False:
                    os.mkdir("biointerpret/feature_files/{}".format(run))

                os.system("cp {}/{}/aggre/feature_aggregation.tab biointerpret/feature_files/{}".format(segwayruns_dir, run, run))
                os.system("cp {}/{}/sigdist/signal_distribution.tab biointerpret/feature_files/{}".format(segwayruns_dir,run, run))

    os.system("cd biointerpret && python apply_samples.py feature_files/ segway_mnemons")
    ls2 = os.listdir("biointerpret/segway_mnemons/")

    for l in ls2:
        if l in ls1:
            os.system("cp biointerpret/segway_mnemons/{}/mnemonics.txt {}/{}".format(l, segwayruns_dir, l))

    # print("cleaning up the feature_files dir...")
    # ls3 = ["biointerpret/feature_files/"+i for i in os.listdir("biointerpret/feature_files") if os.path.isdir("biointerpret/feature_files/"+i)]
    # for l in ls3:
    #     os.system("rm -r {}".format(l))


def chmm_sigdist(chmmruns_dir, original_files_dir):
    ls0 = os.listdir(chmmruns_dir)
    ls1 = os.listdir(original_files_dir)

    for ct in ls1:
        for run in ls0:
            if ct in run:
                ls2 = os.listdir("{}/{}/".format(chmmruns_dir, run))

                if "concat" in run:
                    for ff in ls2:
                        if "dense.bed" in ff:
                            if "rep1" in ff:
                                segbed = "{}/{}/{}".format(chmmruns_dir, run, ff)
                                gd = "{}/{}/rep1.genomedata".format(original_files_dir, ct)
                                outdir = "{}/{}/{}".format(chmmruns_dir, run, 'sigdist_rep1')
                                if os.path.exists(outdir) == False:
                                    os.system(
                                        'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, outdir))

                            elif "rep2" in ff:
                                segbed = "{}/{}/{}".format(chmmruns_dir, run, ff)
                                gd = "{}/{}/rep2.genomedata".format(original_files_dir, ct)
                                outdir = "{}/{}/{}".format(chmmruns_dir, run, 'sigdist_rep2')
                                if os.path.exists(outdir) == False:
                                    os.system(
                                        'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, outdir))
                                        
                elif "rep1psd1" in run:
                    for ff in ls2:
                        if "dense.bed" in ff:
                            segbed = "{}/{}/{}".format(chmmruns_dir, run, ff)
                            gd = "{}/{}/rep1_psdrep1.genomedata".format(original_files_dir, ct)
                            outdir = "{}/{}/{}".format(chmmruns_dir, run, 'sigdist')
                            if os.path.exists(outdir) == False:
                                os.system(
                                    'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, outdir))

                elif "rep1psd2" in run:
                    for ff in ls2:
                        if "dense.bed" in ff:
                            segbed = "{}/{}/{}".format(chmmruns_dir, run, ff)
                            gd = "{}/{}/rep1_psdrep2.genomedata".format(original_files_dir, ct)
                            outdir = "{}/{}/{}".format(chmmruns_dir, run, 'sigdist')
                            if os.path.exists(outdir) == False:
                                os.system(
                                    'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, outdir))

                elif "rep2psd1" in run:
                    for ff in ls2:
                        if "dense.bed" in ff:
                            segbed = "{}/{}/{}".format(chmmruns_dir, run, ff)
                            gd = "{}/{}/rep2_psdrep1.genomedata".format(original_files_dir, ct)
                            outdir = "{}/{}/{}".format(chmmruns_dir, run, 'sigdist')
                            if os.path.exists(outdir) == False:
                                os.system(
                                    'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, outdir))

                elif "rep2psd2" in run:
                    for ff in ls2:
                        if "dense.bed" in ff:
                            segbed = "{}/{}/{}".format(chmmruns_dir, run, ff)
                            gd = "{}/{}/rep2_psdrep2.genomedata".format(original_files_dir, ct)
                            outdir = "{}/{}/{}".format(chmmruns_dir, run, 'sigdist')
                            if os.path.exists(outdir) == False:
                                os.system(
                                    'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, outdir))

                else:
                    for ff in ls2:
                        if "dense.bed" in ff:
                            if "rep1" in ff:
                                segbed = "{}/{}/{}".format(chmmruns_dir, run, ff)
                                gd = "{}/{}/rep1.genomedata".format(original_files_dir, ct)
                                outdir = "{}/{}/{}".format(chmmruns_dir, run, 'sigdist')
                                if os.path.exists(outdir) == False:
                                    os.system(
                                        'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, outdir))

                            elif "rep2" in ff:
                                segbed = "{}/{}/{}".format(chmmruns_dir, run, ff)
                                gd = "{}/{}/rep2.genomedata".format(original_files_dir, ct)
                                outdir = "{}/{}/{}".format(chmmruns_dir, run, 'sigdist')
                                if os.path.exists(outdir) == False:  
                                    os.system(
                                        'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, outdir))


def chmm_aggr(chmmruns_dir, original_files_dir, gtffile="biointerpret/gencode.v29.primary_assembly.annotation_UCSC_names.gtf"):
    ls0 = os.listdir(chmmruns_dir)
    ls1 = os.listdir(original_files_dir)

    for ct in ls1:
        for run in ls0:
            if ct in run:
                ls2 = os.listdir("{}/{}/".format(chmmruns_dir, run))

                if "concat" in run:
                    for ff in ls2:
                        if "dense.bed" in ff:
                            if "rep1" in ff:
                                if os.path.exists("{}/{}/{}".format(chmmruns_dir, run, "aggre_rep1/feature_aggregation.tab")) == False:
                                    os.system('segtools-aggregation --normalize --mode=gene {} {} --outdir={}'.format(
                                    "{}/{}/{}".format(chmmruns_dir, run, ff),
                                    gtffile, "{}/{}/{}".format(chmmruns_dir, run, "aggre_rep1/")))

                            elif "rep2" in ff:
                                if os.path.exists("{}/{}/{}".format(chmmruns_dir, run, "aggre_rep2/feature_aggregation.tab")) == False:
                                    os.system('segtools-aggregation --normalize --mode=gene {} {} --outdir={}'.format(
                                    "{}/{}/{}".format(chmmruns_dir, run, ff),
                                    gtffile, "{}/{}/{}".format(chmmruns_dir, run, "aggre_rep2/")))

                else:
                    for ff in ls2:
                        if "dense.bed" in ff:
                            if os.path.exists("{}/{}/{}".format(chmmruns_dir, run, "aggre/feature_aggregation.tab")) == False:
                                os.system('segtools-aggregation --normalize --mode=gene {} {} --outdir={}'.format(
                                    "{}/{}/{}".format(chmmruns_dir, run, ff),
                                    gtffile, "{}/{}/{}".format(chmmruns_dir, run, "aggre/")))

def chmm_get_mnem(chmmruns_dir):
    """
    for each run:
        mkdir(feature_files/runname)
        cp run/aggre/.tab feature_files/runname/.tab
        cp run/sigdist/.tab feature_files/runname/.tab
    
    run apply.py mnems

    for each run:
        cp mnems/run/mnem.txt run/mnem.txt
    """

    ls1 = os.listdir(chmmruns_dir)

    for run in ls1:
        print("copying files for {}".format(run))
        if "concat" in run:
                
            if os.path.exists("{}/{}/aggre_rep1/feature_aggregation.tab".format(chmmruns_dir, run)) and os.path.exists("{}/{}/sigdist_rep1/signal_distribution.tab".format(chmmruns_dir, run)):
                if os.path.exists("biointerpret/feature_files/{}_rep1".format(run))==False:
                    os.mkdir("biointerpret/feature_files/{}_rep1".format(run))

                os.system("cp {}/{}/aggre_rep1/feature_aggregation.tab biointerpret/feature_files/{}_rep1".format(chmmruns_dir, run, run))
                os.system("cp {}/{}/sigdist_rep1/signal_distribution.tab biointerpret/feature_files/{}_rep1".format(chmmruns_dir, run, run))

            if os.path.exists("{}/{}/aggre_rep2/feature_aggregation.tab".format(chmmruns_dir, run)) and os.path.exists("{}/{}/sigdist_rep2/signal_distribution.tab".format(chmmruns_dir, run)):
                if os.path.exists("biointerpret/feature_files/{}_rep2".format(run))==False:
                    os.mkdir("biointerpret/feature_files/{}_rep2".format(run))

                os.system("cp {}/{}/aggre_rep2/feature_aggregation.tab biointerpret/feature_files/{}_rep2".format(chmmruns_dir, run, run))
                os.system("cp {}/{}/sigdist_rep2/signal_distribution.tab biointerpret/feature_files/{}_rep2".format(chmmruns_dir, run, run))


        else:

            if os.path.exists("{}/{}/aggre/feature_aggregation.tab".format(chmmruns_dir, run)) and os.path.exists("{}/{}/sigdist/signal_distribution.tab".format(chmmruns_dir, run)):
                if os.path.exists("biointerpret/feature_files/{}".format(run))==False:
                    os.mkdir("biointerpret/feature_files/{}".format(run))
 
                os.system("cp {}/{}/aggre/feature_aggregation.tab biointerpret/feature_files/{}".format(chmmruns_dir, run, run))
                os.system("cp {}/{}/sigdist/signal_distribution.tab biointerpret/feature_files/{}".format(chmmruns_dir, run, run))


    os.system("cd biointerpret && python apply_samples.py feature_files/ chmm_mnemons")


    ls2 = os.listdir("biointerpret/chmm_mnemons/")

    for l in ls2:
        if l in ls1 or l.replace("concat_rep1", "concat") in ls1 or l.replace("concat_rep2", "concat") in ls1:
            print("copying files back for {}".format(l))

            if "concat" in l:

                if "rep1" in l:
                    os.system("cp biointerpret/chmm_mnemons/{}/mnemonics.txt {}/{}".format(
                        l, chmmruns_dir,l.replace("_rep1","")))
                    os.system("mv {}/{}/mnemonics.txt {}/{}/mnemonics_rep1.txt".format(
                        chmmruns_dir,l.replace("_rep1",""), chmmruns_dir, l.replace("_rep1","")))
                        
                elif "rep2" in l:
                    os.system("cp biointerpret/chmm_mnemons/{}/mnemonics.txt {}/{}".format(
                        l, chmmruns_dir,l.replace("_rep2","")))
                    os.system("mv {}/{}/mnemonics.txt {}/{}/mnemonics_rep2.txt".format(
                        chmmruns_dir,l.replace("_rep2",""), chmmruns_dir, l.replace("_rep2","")))
            else:
                os.system("cp biointerpret/chmm_mnemons/{}/mnemonics.txt {}/{}".format(l, chmmruns_dir,l))
    
    # print("cleaning up the feature_files dir...")
    # ls3 = ["biointerpret/feature_files/"+i for i in os.listdir("biointerpret/feature_files") if os.path.isdir("biointerpret/feature_files/"+i)]
    # for l in ls3:
    #     os.system("rm -r {}".format(l))


"""
make executable bash using "chmod u+x get_mnemons.sh"

**BASH FILE get_mnemons.sh:
    source ~/miniconda3/etc/profile.d/conda.sh
    activate p27segenv
    python interpret.py sigdist 
    conda deactivate
    conda activate segenv
    python interpret.py feataggr
    python interpret.py get_mnemons
    conda deactivate
"""

def main():
    task = sys.argv[1]
    originalfiles = "files/"
    chmmruns = "chromhmm_runs/"
    segwayruns = "segway_runs/"

    if os.path.exists('biointerpret/gencode.v29.primary_assembly.annotation_UCSC_names.gtf')==False:
        gtf_file()

    if task == "sigdist":
        # try:
        print("running chromhmm signal dist")
        chmm_sigdist(chmmruns_dir=chmmruns, original_files_dir=originalfiles)
            
        # except:
        #     pass

        # try:
        print("running segway signal dist")
        segway_sigdist(segwayruns, originalfiles)
        # except:
        #     pass

    elif task == "feataggr":
        # try:
        print("running segway feature aggregation")
        segway_feataggr(segwayruns)
        # except:
        #     pass

        # try:
        print("running chromhmm feature aggregation")
        chmm_aggr(chmmruns_dir=chmmruns, original_files_dir=originalfiles)
        # except:
        #     pass

    elif task == "mnemon":
        # try:
        print("getting chromhmm mnemons")
        chmm_get_mnem(chmmruns)
        # except:
        #     pass

        # try:
        print("getting segway mnemons")
        segway_get_mnem(segwayruns)
        # except:
        #     pass

if __name__=="__main__":
    main()