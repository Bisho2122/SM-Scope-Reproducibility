import os
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.sparse import coo_matrix
import pickle
import anndata as ad

# Set to your own cluster storage path
CLUSTER_BASE_DIR = os.environ.get("SPATIAL_OMICS_CLUSTER_DIR", "/path/to/your/cluster/storage/Spatial_omics")

def Get_cell_ion_mat_from_coef(coef_pkl,combined_metabo_adata):
    pixel_ion_mat = combined_metabo_adata[coef_pkl['row_idx'].astype('str'),:].X
    pixel_ion_mat = sparse.csr_matrix(pixel_ion_mat)

    cell_pixel_mat = coef_pkl['coef_mat'].T

    cell_ion = cell_pixel_mat.dot(pixel_ion_mat)
    cell_ion = cell_ion.todense()

    return np.array(cell_ion)

def create_multi_layer_sc_metabo_adata(deconv_coef_pkl_dict, cells_transc_adata, combined_metabo_adata,
                                       base_coef_method = 'inverse_sampling_proportion_weighted_by_sampling_specificity_corrected'):

    cell_ion = Get_cell_ion_mat_from_coef(coef_pkl=deconv_coef_pkl_dict[base_coef_method],
                                          combined_metabo_adata=combined_metabo_adata)
    transc_adata_base = cells_transc_adata[deconv_coef_pkl_dict[base_coef_method]['col_idx'],:]

    sc_metabo_adata = ad.AnnData(
        X=cell_ion,
        var=combined_adata.var,
        obs=transc_adata_base.obs
    )
    sc_metabo_adata.obsm['spatial'] = transc_adata_base.obsm['spatial']
    sc_metabo_adata.layers[base_coef_method] = sc_metabo_adata.X

    other_methods = deconv_coef_pkl_dict
    other_methods.pop(base_coef_method)

    for method_name in list(other_methods.keys()):
        cell_ion_mat = Get_cell_ion_mat_from_coef(coef_pkl=deconv_coef_pkl_dict[method_name],
                                          combined_metabo_adata=combined_metabo_adata)

        sc_metabo_adata.layers[method_name] = cell_ion_mat

    return sc_metabo_adata

combined_adata = ad.read_h5ad(f"{CLUSTER_BASE_DIR}/Results/Upstream/SpatialData/M1_S1_S1_final_pixel_all_db_adata.h5ad")

cells_transc_adata = ad.read_h5ad('../../../Data/Figure 1/Resolve_adata.h5ad')
cells_transc_adata.obsm['spatial'] = cells_transc_adata.obs[['x','y']]

deconv_coef_dict = {}
coef_res_path = f"{CLUSTER_BASE_DIR}/Results/Upstream/Pixel_cell_deconvolution/Deconv_Coef_matrices/"
methods_list = ['biggest_overlap_raw','inverse_sampling_proportion_weighted_by_overlap_area',
                'inverse_sampling_proportion_weighted_by_sampling_specificity_corrected',
                'mass_density_weighted_by_overlap_areas',
                'mean_sampling_proportion_corrected',
                'weighted_by_overlap_area',
                'weighted_mean_sampling_area_mark_cell_overlap_int']

for m in methods_list:
    pkl_file = coef_res_path + m + '.pkl'
    with open(pkl_file, 'rb') as f:
        pkl_data = pickle.load(f)
    deconv_coef_dict[m] = pkl_data

final_combined_sc_metabo_adata = create_multi_layer_sc_metabo_adata(deconv_coef_pkl_dict= deconv_coef_dict,
                                                                    cells_transc_adata=cells_transc_adata,
                                                                    combined_metabo_adata=combined_adata,
                                                                    base_coef_method = 'inverse_sampling_proportion_weighted_by_sampling_specificity_corrected')

final_combined_sc_metabo_adata.write_h5ad('../../../Data/Figure 1/M1_S1_S1_cell_ion_adata_w_multiple_deconv.h5ad')