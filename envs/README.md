# Environments

Package versions for the environments this pipeline was actually run in.

- `r_package_versions.csv` / `r_sessionInfo.txt` — R packages used for `Analysis_scripts/Downstream/` and `Fig_scripts/` (R 4.1.2, run locally, not on the cluster).
- Conda environments used for `Analysis_scripts/Upstream/`, one per pipeline stage, all run on the cluster (Linux). Recreate one with `conda create --name <env> --file envs/<file>.yml`:

  | Env file | Pipeline stage |
  |-|-|
  | `env_sm38.yml` | `METASPACE_annotations/` |
  | `env_allensdk_exact.yml` | `Allen_Atlas_to_cells/` (SVG/ontology registration, cell-atlas aggregation) |
  | `env_Image_analysis.yml` | `Allen_Atlas_to_cells/Affine_transform_atlas_RMSE.py` |
  | `env_spatialdata_exact.yml` | `SpatialData/`, `Pixel_Cell_Deconvolution/Sdata_pixel_cell_aggregate.py`, `Build_multimodal_objects/Create_main_METASPACE_adata.py` |
  | `env_spacem.yml` | `Pixel_Cell_Deconvolution/Run_deconv_method.py` |
  | `env_spatialglue_exact.yml` | `Downstream/SpatialGlue/` (Python side) |
