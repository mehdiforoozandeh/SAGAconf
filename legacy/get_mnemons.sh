#!/bin/bash

source  ~/miniconda3/etc/profile.d/conda.sh
conda activate p27segenv
python _interpret.py sigdist 
conda deactivate
conda activate segenv
python _interpret.py feataggr
python _interpret.py mnemon
conda deactivate

