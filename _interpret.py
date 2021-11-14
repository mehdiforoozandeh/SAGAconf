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
    '''
    segwayOutput/exp_name
    segwayOutput/exp_name/feature_aggegation.tab
    segwayOutput/exp_name/signal_distribution.tab
    '''
    os.mkdir('label_interpretation/segwayOutput/exp_name/')
    # os.mkdir('label_interpretation/segwayOutput/exp_name/plots/')

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

def run_apply():
    '''
    run apply_samples.py as follows:
        python apply_samples.py output_folder
    '''
    os.system('python apply_samples.py model_outputs')

def gather_output(exp_name):
    '''
    clean up the results,
    relocate important output files to an isolated location
    '''
    pass