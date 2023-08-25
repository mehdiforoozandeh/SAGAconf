
#!/bin/bash

# Change to the example directory
cd "$(dirname "$0")"
echo "Changed to the example directory"

# Download ChromHMM and unzip it
echo "Downloading ChromHMM..."
curl -o ./ChromHMM.zip http://compbio.mit.edu/ChromHMM/ChromHMM.zip
echo "Unzipping ChromHMM..."
unzip ./ChromHMM.zip

# Generate test annotations from ChromHMM
echo "Generating test annotations from ChromHMM..."
java -mx1600M -jar ./ChromHMM/ChromHMM.jar LearnModel -printposterior ./ChromHMM/SAMPLEDATA_HG18 ./ChromHMM/OUTPUTSAMPLE 10 hg18

# Move posterior files for GM12878_chr11 and K562_chr11 into separate directories for base and verification replicates
echo "Moving posterior files for GM12878_chr11 and K562_chr11 into separate directories for base and verification replicates..."
mkdir -p ./ChromHMM/OUTPUTSAMPLE/base ChromHMM/OUTPUTSAMPLE/verif && \
mv ./ChromHMM/OUTPUTSAMPLE/POSTERIOR/GM12878_10_chr11_posterior.txt ./ChromHMM/OUTPUTSAMPLE/base && \
mv ./ChromHMM/OUTPUTSAMPLE/POSTERIOR/K562_10_chr11_posterior.txt ./ChromHMM/OUTPUTSAMPLE/verif && \
rm -r ./ChromHMM/OUTPUTSAMPLE/POSTERIOR

# Parse posteriors from ChromHMM/OUTPUTSAMPLE/base and ChromHMM/OUTPUTSAMPLE/verif into the standard format required by SAGAconf
echo "Parsing posteriors from ChromHMM/OUTPUTSAMPLE/base and ChromHMM/OUTPUTSAMPLE/verif into the standard format required by SAGAconf..."
python ../SAGAconf_parser.py --saga chmm ./ChromHMM/OUTPUTSAMPLE/base 200 ./ChromHMM/OUTPUTSAMPLE/base
python ../SAGAconf_parser.py --saga chmm ./ChromHMM/OUTPUTSAMPLE/verif 200 ./ChromHMM/OUTPUTSAMPLE/verif

# Run SAGAconf to obtain a full reproducibility report with default parameters
echo "Running SAGAconf to obtain a full reproducibility report with default parameters..."
python ../SAGAconf.py ./ChromHMM/OUTPUTSAMPLE/base/parsed_posterior.bed ./ChromHMM/OUTPUTSAMPLE/verif/parsed_posterior.bed ./sagaconf_base
