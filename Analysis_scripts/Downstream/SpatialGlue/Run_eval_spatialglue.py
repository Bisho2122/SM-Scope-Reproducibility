import scib
import scanpy as sc
import pandas as pd
import pickle
import os
import argparse

# Initialize parser
parser = argparse.ArgumentParser(description="SpatialGlue param input")

# Add arguments
parser.add_argument("--input", type=str, help="Spatial glue results")

# Parse arguments
args = parser.parse_args()

def Run_scib_per_group(adata_obj, label_grp):
    # Ensure label_grp is in adata.obs
    if label_grp not in adata_obj.obs:
        raise ValueError(f"Column '{label_grp}' not found in adata.obs")

    # Remove NaN values in the label column
    adata_scib = adata_obj[~adata_obj.obs[label_grp].isna()].copy()

    adata_scib.obs[label_grp] = adata_scib.obs[label_grp].astype("category")

    # Validate required keys
    if "SpatialGlue" not in adata_scib.obsm:
        raise ValueError("Embedding 'SpatialGlue' not found in adata.obsm")

    if "SG_leiden" not in adata_scib.obs:
        raise ValueError("Column 'SG_leiden' not found in adata.obs")

    # Compute metrics
    asw = scib.metrics.silhouette(adata_scib, label_key=label_grp, embed='SpatialGlue')
    clisi = scib.metrics.clisi_graph(adata_scib, label_key=label_grp,
                                     use_rep='SpatialGlue', type_='embed')
    ari = scib.metrics.ari(adata_scib, cluster_key="SG_leiden", label_key=label_grp)
    nmi = scib.metrics.nmi(adata_scib, cluster_key="SG_leiden", label_key=label_grp)

    # Store results in a DataFrame
    final_df = pd.DataFrame({
        "Label": [label_grp] * 4,
        "Metric": ["asw", "clisi", "ari", "nmi"],
        "Value": [asw, clisi, ari, nmi]
    })

    return final_df

def Run_scib_multiple_groups(adata_obj, label_grps):
    results = []  # List to store DataFrames

    for label_grp in label_grps:
        try:
            df = Run_scib_per_group(adata_obj, label_grp)  # Run for each label group
            results.append(df)
        except ValueError as e:
            print(f"Skipping '{label_grp}': {e}")  # Handle errors gracefully

    # Combine all results into a single DataFrame
    final_df = pd.concat(results, ignore_index=True) if results else pd.DataFrame()

    return final_df

#Main
adata_metabo = sc.read_h5ad('../../../Data/Figure 1/M1_S1_S1_cell_ion_adata_w_multiple_deconv.h5ad')

with open(args.input, "rb") as f:
    sg_res = pickle.load(f)

grp_metadata = pd.read_csv("../../../Data/Figure 6/Cell_grps_metadata.csv", index_col = 0)
adata_metabo.obs = adata_metabo.obs.join(grp_metadata, how="left")

#Adding metadata and spatialglue results
adata = adata_metabo.copy()
adata.obsm['emb_latent_omics1'] = sg_res['emb_latent_omics1'].copy()
adata.obsm['emb_latent_omics2'] = sg_res['emb_latent_omics2'].copy()
adata.obsm['SpatialGlue'] = sg_res['SpatialGlue'].copy()
adata.obsm['alpha'] = sg_res['alpha']
adata.obsm['alpha_omics1'] = sg_res['alpha_omics1']
adata.obsm['alpha_omics2'] = sg_res['alpha_omics2']

adata.obs['alpha_RNA'] = sg_res['alpha'][:,0]
adata.obs['alpha_metabo'] = sg_res['alpha'][:,1]
adata.obs['alpha_metabo_spatial'] = sg_res['alpha_omics2'][:,0]
adata.obs['alpha_metabo_feat'] = sg_res['alpha_omics2'][:,1]
adata.obs['alpha_RNA_spatial'] = sg_res['alpha_omics1'][:,0]
adata.obs['alpha_RNA_feat'] = sg_res['alpha_omics1'][:,1]

adata.obsm['spatial'] = adata.obsm['spatial'].to_numpy()

#Clustering using optimal res
adata_clust = adata[~adata.obs["Region"].isna()].copy()
optim_resol = scib.metrics.cluster_optimal_resolution(adata_clust,
                                                      label_key='Region',
                                                      use_rep='SpatialGlue',
                                                      cluster_key="SG_leiden",
                                                      return_all = True)


sc.pp.neighbors(adata, use_rep='SpatialGlue', n_neighbors=15)
sc.tl.leiden(adata, key_added='SG_leiden', resolution= optim_resol[0])

#scIB pipeline

all_cell_grps = ["Cell_type_coarse", "Region", "CT_Region"]

metrics_results = Run_scib_multiple_groups(adata_obj=adata, label_grps=all_cell_grps)

#Save output
file_name = os.path.splitext(os.path.basename(args.input))[0]

# Raw per-parameter-combination SpatialGlue eval output from the SLURM sweep -- not checked
# into Data/, matches root_dir in SpatialGlue_param_optim.R. Set to your own cluster storage path.
CLUSTER_BASE_DIR = os.environ.get("SPATIAL_OMICS_CLUSTER_DIR", "/path/to/your/cluster/storage/Spatial_omics")
out_metrics_path = f"{CLUSTER_BASE_DIR}/Results/Downstream/SpatialGlue/SpatialGlue_eval/SCIB_metrics/{file_name}.csv"
out_meta_path = f"{CLUSTER_BASE_DIR}/Results/Downstream/SpatialGlue/SpatialGlue_eval/SG_clusts_metadata/{file_name}.csv"
out_loss_path = f"{CLUSTER_BASE_DIR}/Results/Downstream/SpatialGlue/SpatialGlue_eval/SG_loss/{file_name}.csv"

metrics_results.to_csv(out_metrics_path)
adata.obs.to_csv(out_meta_path)
pd.DataFrame(sg_res['loss_dict']).to_csv(out_loss_path)