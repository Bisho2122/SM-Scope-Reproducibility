# Data

This repo's data is split across **two Zenodo records**:

- **Record 1 (this folder, `Data/`)**: processed/intermediate objects (Seurat objects, AnnData files, spatial niche/region assignments, MISTY/DIALOGUE/SpatialGlue results, KEGG/Reactome/RHEA reaction-network objects). This is everything `Analysis_scripts/Downstream/` and `Fig_scripts/` need; **required** to reproduce the paper's analyses and figures.
- **Record 2 (raw + upstream cluster intermediates)**: raw instrument data and the `Analysis_scripts/Upstream/` pipeline's intermediate artifacts. Only needed if you want to regenerate Record 1's objects from scratch; not required to reproduce the paper's downstream analyses/figures.

`Data/` is **not tracked in git** (see `.gitignore`).

## Record 1: getting `Data/`

**DOI: TODO**
**Record ID: TODO**

Once the record exists, download and unpack it from the repo root with:

```bash
zenodo_download_data.sh <zenodo_record_id>
```

This fetches each archive from the Zenodo record, verifies its checksum against the one Zenodo reports, and unpacks it into this folder, reproducing the layout below. Requires `curl`, `jq`, and `unzip`.

### Expected layout

```text
Data/
├── Figure 1/    # inputs/outputs for Figure 1; also the primary input for most downstream scripts
├── Figure 2/
├── Figure 3/
├── Figure 4/
├── Figure 5/
├── Figure 6/
├── General/     # cross-figure reference data (ontology tables, ion/molecule name mapping)
└── *.rds / *.RData
```

## Record 2: getting the raw/upstream data

Only needed if you're re-running `Analysis_scripts/Upstream/` scripts or the SpatialGlue parameter sweep from scratch.

**DOI: TODO**
**Record ID: TODO**

Download it into a local directory of your choice, then point `Upstream/`/sweep scripts at it via one environment variable:

```bash
zenodo_download_data.sh <zenodo_record_id> /path/to/local/upstream_data
export SPATIAL_OMICS_CLUSTER_DIR=/path/to/local/upstream_data
```

Every script that reads this data falls back to a `/path/to/your/cluster/storage/Spatial_omics`-style placeholder if the env var isn't set, so this is the only setup step needed before running them.

### Expected layout

```text
<SPATIAL_OMICS_CLUSTER_DIR>/
├── Data/
│   ├── SpatialData_input/
│   │   ├── 32753-Slide1_A1_DAPI.tiff
│   │   ├── M1_Sand1_Sec1_cropped.tif
│   │   ├── M1_Sand1_sec1_ion.tif
│   │   ├── Section_68_annotated.jpg
│   │   ├── M1_Sand1_Sec1_to_ISS.txt
│   │   ├── Ion_to_M1_Sand1_Sec1.txt
│   │   ├── Final_section_68_whole_to_cells_HPF_mat.txt
│   │   └── updated_section_68_annots.geojson
│   └── Data_Delivery_32753_S1_D1/32753-Slide1_submission/
│       └── 32753_Slide1_A1_spots.csv
└── Results/
    ├── QuPath/
    │   └── 32753-Slide1_A1_DAPI-detections.tif
    ├── Baysor_segmentation/Baysor_run_w_segmask_out_dev/
    │   └── segmentation.csv   # (plus other Baysor output files, kept for provenance but not read by the current pipeline)
    ├── Upstream/
    │   ├── SpatialData/
    │   │   ├── Final_SpatialData_obj_April_22_2025/    # a SpatialData store (directory), not a single file -- extract fully before use
    │   │   ├── M1_S1_S1_final_pixel_all_db_adata.h5ad
    │   │   ├── updated_geopandas_intersect_cells_atlas.csv
    │   │   └── updated_geopandas_intersect_cells_pixels.csv
    │   ├── Pixel_cell_deconvolution/Deconv_Coef_matrices/
    │   │   └── <method-name>.pkl files
    │   ├── Combined_obj_for_integration/
    │   │   ├── M1_S1_S1_combined_adata_var_df.csv
    │   │   └── M1_S1_S1_combined_adata_var_df_w_new_spectral.csv
    │   └── METASPACE_annotations/
    │       ├── 2022-01-11_23h33m36s.csv
    │       └── Annotation_QC/
    └── Downstream/SpatialGlue/
        ├── SpatialGlue_results/
        │   └── <per-param-combo>.pkl (many, one per sweep parameter combination)
        ├── updated_dim50_k15_metriceuclidean_wf5511.pkl   # canonical best-param model, kept separate from the sweep
        └── SpatialGlue_eval/{SCIB_metrics,SG_clusts_metadata,SG_loss}/*.csv
```
