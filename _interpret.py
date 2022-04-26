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

def get_mnem():
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

def run_sigdist_concat():
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