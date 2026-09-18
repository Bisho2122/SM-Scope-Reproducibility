#!/bin/bash

# Set these to your own local paths (see README.md for the expected repo layout)
R_LIBRARY_DIR="${SPATIAL_OMICS_R_LIBRARY_DIR:-/path/to/your/R/library}"
REPO_DIR="${SPATIAL_OMICS_REPO_DIR:-/path/to/your/SM-Scope-Reproducibility}"
DIALOGUE_SIF="${DIALOGUE_SINGULARITY_IMAGE:-/path/to/your/DIALOGUE_singularity/DIALOGUE_image.sif}"

singularity exec \
  --bind "$R_LIBRARY_DIR":/mnt/library \
  --bind "$REPO_DIR/Analysis_scripts/Downstream/DIALOGUE":/mnt/scripts \
  --bind "$REPO_DIR/Data/Figure 5":/mnt/inputs \
  --bind "$REPO_DIR/Data/Figure 5":/mnt/outputs \
  "$DIALOGUE_SIF" \
  Rscript /mnt/scripts/Run_DIALOGUE_singularity.R /mnt/inputs/DIALOGUE_input.rds /mnt/outputs/updated_DIALOGUE_res
