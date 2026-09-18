#!/bin/bash

#SBATCH -J Run_SG
#SBATCH --time=01-00:00:00
#SBATCH --mem=30000M
#SBATCH -n 2
#SBATCH -c 2
#SBATCH --mail-type=FAIL,END

dim="$1"
k="$2"
metric="$3"
wf="$4"

# Run the SpatialGlue python script
python Run_spatialglue.py --dim_out "$dim" --k_feat "$k" --metric_feat "$metric" --wf "$wf" --epochs 600