#!/bin/bash

conda activate p27segenv
python _interpret.py sigdist 
conda deactivate
conda activate segenv
python _interpret.py feataggr
python _interprete.py mnemon
conda deactivate

