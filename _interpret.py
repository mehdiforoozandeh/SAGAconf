import wget, os

def gtf_file():
    if not os.path.isfile('label_interpretation/gencode.v29.primary_assembly.annotation_UCSC_names.gtf'):
        url='https://www.encodeproject.org/files/\
        gencode.v29.primary_assembly.annotation_UCSC_names/\
            @@download/gencode.v29.primary_assembly.annotation_UCSC_names.gtf.gz'

        filename = wget.download(url, out='label_interpretation')
        os.system('gzip -d {}'.format(filename))

        return filename.replace('.gz', '')
    
def feature_aggreg(exp_name, segbed, gtf):
    '''
    use segtools to run segtools feat_aggreg
    '''
    outdir = 'label_interpretation/segwayOutput/{}'.format(exp_name)
    
    os.system('segtools-aggregation --normalize {} {} --outdir={}'.format(segbed, gtf, outdir))

def signal_dist(exp_name, segbed, gd):
    '''
    use segtools to run segtools signal_distribution
    '''
    outdir = 'label_interpretation/segwayOutput/{}'.format(exp_name)

    os.system('segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, outdir))

def create_input_dir():
    '''
    segwayOutput/exp_name
    segwayOutput/exp_name/feature_aggegation.tab
    segwayOutput/exp_name/signal_distribution.tab
    segwayOutput/exp_name/plots/
    '''
    pass

def trackname_guide():
    '''
    generate 
    segwayOutput/trackname_assays.txt
    '''

def run_apply(exp_name):
    '''
    run apply_samples.py as follows:
        python apply_samples.py output_folder
    '''
    pass

def run_plots(exp_name):
    pass

def gather_output(exp_name):
    '''
    clean up the results,
    relocate important output files to an isolated location
    '''
    pass