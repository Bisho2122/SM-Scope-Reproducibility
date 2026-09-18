#!/bin/bash

#SBATCH -J Aggregate_atlas_sdata
#SBATCH --time=02-00:00:00
#SBATCH --mem=150000M
#SBATCH -n 2
#SBATCH -c 2
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH -e slurm-%x.err
#SBATCH -o slurm-%x.out

python atlas_aggregate_job.py
