import os
import torch
import pandas as pd
import scanpy as sc

import argparse
import pickle

import SpatialGlue.preprocess as sg_pp
import SpatialGlue.SpatialGlue_pyG as sg_model

# Initialize parser
parser = argparse.ArgumentParser(description="SpatialGlue param input")

# Add arguments
parser.add_argument("--dim_out", type=int, default=50, help="Spatial Glue embedding out dimension")
parser.add_argument("--n_neigh", type=int, default=3, help="n neighbors for spatial graph")
parser.add_argument("--k_feat", type=int, default=20, help="k neighbors for feature graph")
parser.add_argument("--metric_feat", type=str, default="correlation", help="KNN metric for feature graph")
parser.add_argument("--pca_comps", type=int, default=30, help="Number of PCA components for preprocessing")
parser.add_argument("--epochs", type=int, default=600, help="Number of epochs for training")
parser.add_argument("--wf", type=str, default="1511", help="String of weight factors for 4 losses in spatial glue")

# Parse arguments
args = parser.parse_args()

# Read data
adata_RNA = sc.read_h5ad('../../../Data/Figure 1/Resolve_adata.h5ad')
adata_metabo = sc.read_h5ad('../../../Data/Figure 1/M1_S1_S1_cell_ion_adata_w_multiple_deconv.h5ad')

#config
data_type = 'RNA_metabo'
random_seed = 42
sg_pp.fix_seed(random_seed)

dim_out = args.dim_out
n_neigh = args.n_neigh
k_feat = args.k_feat
metric_feat = args.metric_feat
pca_comps = args.pca_comps
n_epochs = args.epochs
loss_weight_factors = [int(i) for i in list(args.wf)]

#Preprocessing

common_cells = adata_metabo.obs_names.intersection(adata_RNA.obs_names)
adata_metabo = adata_metabo[common_cells,]
adata_RNA = adata_RNA[common_cells,]

sc.pp.calculate_qc_metrics(adata_metabo, inplace=True, percent_top=(10,50,100))
sc.pp.filter_cells(adata_metabo, min_counts=600)
sc.pp.filter_cells(adata_metabo, max_counts=10000)
sc.pp.filter_genes(adata_metabo, min_cells=10)
sc.pp.normalize_total(adata=adata_metabo,target_sum=1e4)

adata_omics_2 = adata_metabo.copy()
adata_omics_2.raw = adata_metabo
adata_omics_2.X = adata_omics_2.layers['inverse_sampling_proportion_weighted_by_sampling_specificity_corrected']
sc.pp.normalize_total(adata=adata_omics_2, target_sum=1e4)
sc.pp.log1p(adata_omics_2)
sc.pp.scale(adata_omics_2)

adata_omics_2.obsm['feat'] = sg_pp.pca(adata_omics_2, n_comps=pca_comps)

#preprocess RNA

sc.pp.filter_genes(adata_RNA, min_cells=10)
sc.pp.normalize_total(adata_RNA, target_sum=1e4)
sc.pp.log1p(adata_RNA)
sc.pp.scale(adata_RNA)

adata_RNA.obsm['feat'] = sg_pp.pca(adata_RNA, n_comps=pca_comps)
adata_omics_1 = adata_RNA.copy()

adata_omics_1.obsm['spatial'] = adata_omics_2.obsm['spatial']

#Construct neighbor graphs

sg_data = sg_pp.construct_neighbor_graph(adata_omics1=adata_omics_1,
                                         adata_omics2=adata_omics_2,
                                         n_neighbors=n_neigh,
                                         datatype=data_type,
                                         k = k_feat,
                                         metric = metric_feat)


#Model training
torch_device = torch.device('cuda:2' if torch.cuda.is_available() else 'cpu')

model = sg_model.Train_SpatialGlue(sg_data, datatype=data_type,
                                   device=torch_device,
                                   random_seed= random_seed,
                                   epochs=n_epochs,
                                   dim_output=dim_out,
                                   weight_factors=loss_weight_factors,
                                   dim_input=50)

sg_output = model.train()

out_filename = f"dim{args.dim_out}_k{args.k_feat}_metric{args.metric_feat}_wf{args.wf}.pkl"

# Raw per-parameter-combination SpatialGlue output from the SLURM sweep (SpatialGlue_job_parallel.sh)
# -- not checked into Data/, matches root_dir in SpatialGlue_param_optim.R. Set to your own cluster storage path.
CLUSTER_BASE_DIR = os.environ.get("SPATIAL_OMICS_CLUSTER_DIR", "/path/to/your/cluster/storage/Spatial_omics")
out_path = f"{CLUSTER_BASE_DIR}/Results/Downstream/SpatialGlue/SpatialGlue_results/{out_filename}"
os.makedirs(os.path.dirname(out_path), exist_ok=True)

# Save to file
with open(out_path, "wb") as f:
    pickle.dump(sg_output, f)

print("All Done")
