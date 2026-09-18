import os
from metaspace_converter import metaspace_to_anndata
from metaspace import SMInstance
import anndata as ad
import numpy as np
import pandas as pd

# Set to your own cluster storage path (matches Create_multi_layer_RNA_metabo_adata.py's read of the same object)
CLUSTER_BASE_DIR = os.environ.get("SPATIAL_OMICS_CLUSTER_DIR", "/path/to/your/cluster/storage/Spatial_omics")

def combine_anndata_list(adata_list):
    # Check that all AnnData objects have the same observations
    obs_names = adata_list[0].obs_names
    for adata in adata_list:
        assert (adata.obs_names == obs_names).all(), "All AnnData objects must have the same observations!"

    # Get the union of all variable names
    all_vars = [i.var_names for i in adata_list]
    all_vars = [item for sublist in all_vars for item in sublist]
    all_vars = pd.Index(all_vars)

    # Initialize the combined data matrix
    combined_data = np.zeros((adata_list[0].shape[0], len(all_vars)))

    # Fill the matrix with data from each AnnData object
    for adata in adata_list:
        adata_indices = all_vars.get_indexer(adata.var_names)
        combined_data[:, adata_indices] += adata.X.toarray() if isinstance(adata.X, np.matrix) else adata.X

    # Create the combined AnnData object

    var_meta = pd.concat([adata.var for adata in adata_list])
    combined_adata = ad.AnnData(
        X=combined_data,
        obs=adata_list[0].obs.copy(),  # Use obs from the first AnnData
        var=var_meta  # Union of variables

    )
    combined_adata.uns = adata_list[0].uns.copy()
    combined_adata.obsm = adata_list[0].obsm.copy()

    return combined_adata

sm = SMInstance()
M1_S1_S1_ds_id= '2022-01-11_23h33m36s'
ds = sm.dataset(id = M1_S1_S1_ds_id)
databases_to_consider = ds.database_details

annot_QC = []
for db in ['CoreMetabolome','HMDB','CHEBI','RAMP_DB', 'Refmet']:
    annot_qc_path = f"{CLUSTER_BASE_DIR}/Results/Upstream/METASPACE_annotations/Annotation_QC/Annot_QC_2022-01-11_23h33m36s_{db}.csv"
    qc_data = pd.read_csv(annot_qc_path)
    annot_QC.append(qc_data)

all_annot_QC = pd.concat(annot_QC)
all_annot_QC = all_annot_QC[all_annot_QC.Sparsity < 0.99]
all_annot_QC = all_annot_QC[all_annot_QC.off_sample_label == 'on']
all_annot_QC = all_annot_QC.drop_duplicates(['formula','adduct'])
all_annot_QC['ion'] = all_annot_QC['formula'] + all_annot_QC['adduct']

filtered_adata_list = []
for db in databases_to_consider:
    qc_db = all_annot_QC[all_annot_QC.DB_name == db.name]
    adata_obj = metaspace_to_anndata(dataset_id=M1_S1_S1_ds_id, fdr = 1,database=(db.name,db.version),
                                    sm=sm)
    adata_obj = adata_obj[:,qc_db.ion]
    filtered_adata_list.append(adata_obj)

adata_list_pass_QC = [i[:,i.var_names.isin(list(all_annot_QC.ion))] for i in filtered_adata_list]

unique_var_names = []
unique_adata_list = []
for i, adata in enumerate(adata_list_pass_QC):
    new_adata = adata.copy()
    current_vars = set(adata.var_names)
    if len(unique_var_names) == 0:
        unique_var_names.extend(list(current_vars))
    else:
        overlapping_vars = current_vars.intersection(set(unique_var_names))
        unique_vars = current_vars - overlapping_vars
        new_adata = new_adata[:,new_adata.var_names.isin(list(unique_vars))]
        unique_var_names.extend(list(current_vars))

    unique_adata_list.append(new_adata)

final_combined_adata = combine_anndata_list(unique_adata_list)
final_combined_adata.var = pd.merge(left=final_combined_adata.var, right = all_annot_QC,
                                    how='left', left_on=['formula','adduct'],right_on=['formula','adduct'])
final_combined_adata.var = final_combined_adata.var.set_index('ion_y')

final_combined_adata.write_h5ad(f"{CLUSTER_BASE_DIR}/Results/Upstream/SpatialData/M1_S1_S1_final_pixel_all_db_adata.h5ad")