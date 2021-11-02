## segway 3.0.3 run ee9e72ba38f511ec96cc509a4c776cf2 at 2021-10-29 13:22:26.417295

cd "/scratch/mforooz/segwayconf/GM12878/SAGAconf/segway_runs"
"/home/mforooz/miniconda3/envs/segenv/bin/segway" "train" "--include-coords=encodePilotRegions.hg19.bed" "--num-instances=10" "--track-weight=0.0008333333333333334" "--segtransition-weight-scale=0.01" "--ruler-scale=1000" "--prior-strength=1" "--resolution=1000" "--num-labels=14" "rep2.genomedata" "train_rep2_14labels_res1000"
