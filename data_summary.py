import pandas as pd
import sys, os

def get_summary(download_dir):
    summary = []
    celltypes = [ct for ct in os.listdir(download_dir) if os.path.isdir(download_dir+"/"+ct)]
    for ct in celltypes:
        tracks = [tr for tr in os.listdir(download_dir+"/"+ct) if os.path.isdir(download_dir+"/"+ct+"/"+tr)]
        for tr in tracks:
            tfmd = pd.read_csv(download_dir+"/"+ct+'/'+tr+'/track_files_metadata.csv')
            tfmd.index = list(tfmd['Unnamed: 0'])
            tfmd = tfmd.drop('Unnamed: 0', axis=1)
            summary.append(
                [
                    ct, tfmd.loc["assay", "rep1_alig"], 
                    str(", ".join([
                        tfmd.loc["biosample", "rep1_alig"]+" (rep. 1)", 
                        tfmd.loc["biosample", "rep2_alig"]+" (rep. 2)"
                        ])),
                    str(", ".join([
                        tfmd.loc["accession", "rep1_alig"], tfmd.loc["accession", "rep1_fcoc"], 
                        tfmd.loc["accession", "rep2_alig"], tfmd.loc["accession", "rep2_fcoc"]
                        ]))
                ]
            )
    summary = pd.DataFrame(summary, columns=['Celltype', 'experiment', 'biosample (isoreplicate)', 'files'])
    return summary

if __name__ == "__main__":
    download_dir = sys.argv[1]
    summary = get_summary(download_dir)
    summary.to_csv("data_summary.csv")
