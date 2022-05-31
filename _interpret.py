import requests, os

def gtf_file(gtf_filename = 'label_interpretation/gencode.v29.primary_assembly.annotation_UCSC_names.gtf'): 
    if not os.path.isfile(gtf_filename):
        url='https://www.encodeproject.org/files/\
        gencode.v29.primary_assembly.annotation_UCSC_names/\
            @@download/gencode.v29.primary_assembly.annotation_UCSC_names.gtf.gz'

        r = requests.get(url, allow_redirects=True)
        open(gtf_filename + '.gz', 'wb').write(r.content)

        os.system('gzip -d {}'.format(gtf_filename+ '.gz'))

        return gtf_filename

def make_temp_bed_for_concat(segbed_file):
    "for running segtools on concat runs, we need to convert the bed files to sth with normal chr names"
    write_obj = open(segbed_file.replace(".bed", "_temp.bed"), 'w')
    with open(segbed_file, 'r') as sfn:
        lines = sfn.readlines()
        for l in lines:
            if l[:3] == "chr":
                l = l.split('\t')
                l[0] = l[0][:-2]
                write_obj.write("\t".join(l))
            else:
                write_obj.write(l)
    write_obj.close()

def create_input_dir(exp_name):
    ''' MUST CONTAIN:
    segwayOutput/exp_name
    segwayOutput/exp_name/feature_aggegation.tab
    segwayOutput/exp_name/signal_distribution.tab
    segwayOutput/exp_name/segway.bed
    segwayOutput/exp_name/{gtf_file}
    segwayOutput/exp_name/{genomedata_file}
    '''
    os.mkdir('label_interpretation/segwayOutput/exp_name/')

def feature_aggreg(exp_name, segbed, gtf):
    '''
    use segtools to run segtools feat_aggreg
    '''
    outdir = 'label_interpretation/segwayOutput/{}'.format(exp_name)

    os.system('segtools-aggregation --normalize --mode=gene {} {} --outdir={}'.format(segbed, gtf, outdir))

def signal_dist(exp_name, segbed, gd):
    '''
    use segtools to run segtools signal_distribution
    '''
    outdir = 'label_interpretation/segwayOutput/{}'.format(exp_name)

    os.system('segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, outdir))

def segway_get_mnem():
    """
    for each run:
        mkdir(segwayoutput/runname)
        cp run/aggre/.tab segwayoutput/runname/.tab
        cp run/sigdist/.tab segwayoutput/runname/.tab
    
    run apply.py mnems

    for each run:
        cp mnems/run/mnem.txt run/mnem.txt
        """
    ls1 = os.listdir("segway_runs/")
    for run in ls1:
        with open("segway_runs/{}/signal_dist/signal_distribution.tab".format(run), 'r') as file:
            lines = file.readlines()
            lines = "".join(lines)
        if "nan" not in lines:
            if os.path.exists("biointerpret/segwayOutput/{}".format(run))==False:
                os.mkdir("biointerpret/segwayOutput/{}".format(run))
            os.system("cp segway_runs/{}/aggregations/feature_aggregation.tab biointerpret/segwayOutput/{}".format(run, run))
            os.system("cp segway_runs/{}/signal_dist/signal_distribution.tab biointerpret/segwayOutput/{}".format(run, run))
    os.system("cd biointerpret && python apply_samples.py segway_mnemons")
    ls2 = os.listdir("biointerpret/segway_mnemons/classification")
    for l in ls2:
        os.system("cp biointerpret/segway_mnemons/classification/{}/mnemonics.txt segway_runs/{}".format(l,l))

def segway_run_sigdist_concat():
    ls0 = os.listdir("files/")
    ls1 = os.listdir("segway_runs/")
    for i in ls0:
        for j in ls1:
            print(j)
            segbed = "segway_runs/"+j+'/segway.bed'
            if i in j:
                if "concat" in j:
                    if "rep1" in j:
                        gd = "files/"+i+"/concat_rep1.genomedata"
                    elif "rep2" in j:
                        gd = "files/"+i+"/concat_rep2.genomedata"

                    os.system(
                        'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, "segway_runs/"+j+'/signal_dist'))

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
                                os.system(
                                    'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, outdir))

                            elif "rep2" in ff:
                                segbed = "{}/{}/{}".format(chmmruns_dir, run, ff)
                                gd = "{}/{}/rep2.genomedata".format(original_files_dir, ct)
                                outdir = "{}/{}/{}".format(chmmruns_dir, run, 'sigdist_rep2')
                                os.system(
                                    'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, outdir))

                # else:
                #     for ff in ls2:
                #         if "dense.bed" in ff:
                #             if "rep1" in ff:
                #                 segbed = "{}/{}/{}".format(chmmruns_dir, run, ff)
                #                 gd = "{}/{}/rep1.genomedata".format(original_files_dir, ct)
                #                 outdir = "{}/{}/{}".format(chmmruns_dir, run, 'sigdist')
                #                 os.system(
                #                     'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, outdir))

                #             elif "rep2" in ff:
                #                 segbed = "{}/{}/{}".format(chmmruns_dir, run, ff)
                #                 gd = "{}/{}/rep2.genomedata".format(original_files_dir, ct)
                #                 outdir = "{}/{}/{}".format(chmmruns_dir, run, 'sigdist')
                #                 os.system(
                #                     'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, outdir))


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
                                os.system('segtools-aggregation --normalize --mode=gene {} {} --outdir={}'.format(
                                "{}/{}/{}".format(chmmruns_dir, run, ff),
                                 gtffile, "{}/{}/{}".format(chmmruns_dir, run, "aggre_rep1/")))

                            elif "rep2" in ff:
                                os.system('segtools-aggregation --normalize --mode=gene {} {} --outdir={}'.format(
                                "{}/{}/{}".format(chmmruns_dir, run, ff),
                                 gtffile, "{}/{}/{}".format(chmmruns_dir, run, "aggre_rep2/")))

                else:
                    for ff in ls2:
                        if "dense.bed" in ff:
                            os.system('segtools-aggregation --normalize --mode=gene {} {} --outdir={}'.format(
                                "{}/{}/{}".format(chmmruns_dir, run, ff),
                                 gtffile, "{}/{}/{}".format(chmmruns_dir, run, "aggre/")))

def chmm_get_mnem():
    """
    for each run:
        mkdir(segwayoutput/runname)
        cp run/aggre/.tab segwayoutput/runname/.tab
        cp run/sigdist/.tab segwayoutput/runname/.tab
    
    run apply.py mnems

    for each run:
        cp mnems/run/mnem.txt run/mnem.txt
    """

    ls1 = os.listdir("segway_runs/")
    for run in ls1:
        with open("segway_runs/{}/signal_dist/signal_distribution.tab".format(run), 'r') as file:
            lines = file.readlines()
            lines = "".join(lines)
        if "nan" not in lines:
            if os.path.exists("biointerpret/segwayOutput/{}".format(run))==False:
                os.mkdir("biointerpret/segwayOutput/{}".format(run))
            os.system("cp segway_runs/{}/aggregations/feature_aggregation.tab biointerpret/segwayOutput/{}".format(run, run))
            os.system("cp segway_runs/{}/signal_dist/signal_distribution.tab biointerpret/segwayOutput/{}".format(run, run))
    os.system("cd biointerpret && python apply_samples.py segway_mnemons")
    ls2 = os.listdir("biointerpret/segway_mnemons/classification")
    for l in ls2:
        os.system("cp biointerpret/segway_mnemons/classification/{}/mnemonics.txt segway_runs/{}".format(l,l))

def segway_run_sigdist_concat():
    ls0 = os.listdir("files/")
    ls1 = os.listdir("segway_runs/")
    for i in ls0:
        for j in ls1:
            print(j)
            segbed = "segway_runs/"+j+'/segway.bed'
            if i in j:
                if "concat" in j:
                    if "rep1" in j:
                        gd = "files/"+i+"/concat_rep1.genomedata"
                    elif "rep2" in j:
                        gd = "files/"+i+"/concat_rep2.genomedata"

                    os.system(
                        'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, "segway_runs/"+j+'/signal_dist'))
