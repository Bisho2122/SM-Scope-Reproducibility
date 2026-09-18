from SpaceM_Deconv_methods import *
import os
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.sparse import coo_matrix
import xarray as xr
import sys
import pickle

if len(sys.argv) != 2:
    print("Usage: python process_string.py <string>")
    sys.exit(1)

deconv_method = sys.argv[1]

# Set to your own cluster storage path (matches Create_multi_layer_RNA_metabo_adata.py's read of the same intermediates,
# and Sdata_pixel_cell_aggregate.py's write of the geopandas intersect csv)
CLUSTER_BASE_DIR = os.environ.get("SPATIAL_OMICS_CLUSTER_DIR", "/path/to/your/cluster/storage/Spatial_omics")

out_path = f"{CLUSTER_BASE_DIR}/Results/Upstream/Pixel_cell_deconvolution/Deconv_Coef_matrices/"

gpd_pixel_res = pd.read_csv(f"{CLUSTER_BASE_DIR}/Results/Upstream/SpatialData/updated_geopandas_intersect_cells_pixels.csv")

def pivot_wider_and_sparse_mat(df, rows_colname, cols_colname, value_colname):
    rows = df[rows_colname].astype('category').cat.codes
    cols = df[cols_colname].astype('category').cat.codes
    values = df[value_colname]

    row_labels = df[rows_colname].astype('category').cat.categories
    col_labels = df[cols_colname].astype('category').cat.categories

    csr_mat = coo_matrix(
    (values, (rows, cols)),
    shape=(len(df[rows_colname].unique()), len(df[cols_colname].unique()))
    ).tocsr()

    final = {'row_idx' : row_labels,
             'col_idx' : col_labels,
             'mat' : csr_mat}
    return final

overlap_matrix = pivot_wider_and_sparse_mat(df = gpd_pixel_res,
                                            rows_colname='instance',
                                            cols_colname='__index',
                                            value_colname='intersect_areas')

overlap_matrix['mat'] = xr.DataArray(overlap_matrix['mat'].todense(), coords=[overlap_matrix['row_idx'], overlap_matrix['col_idx']])

cell_areas = gpd_pixel_res[['__index','cell_areas']].drop_duplicates()
cell_areas.index = cell_areas['__index']
cell_areas = cell_areas.reindex(overlap_matrix['col_idx'])
assert cell_areas.index.equals(overlap_matrix['col_idx'])
cell_areas = cell_areas['cell_areas'].to_numpy()

pixel_areas = gpd_pixel_res[['instance','instance_areas']].drop_duplicates()
pixel_areas.index = pixel_areas['instance']
pixel_areas = pixel_areas.reindex(overlap_matrix['row_idx'])
assert pixel_areas.index.equals(overlap_matrix['row_idx'])
pixel_areas = pixel_areas['instance_areas'].to_numpy()

# pixel_areas = np.reshape(pixel_areas, (-1,1))

deconv_function = DECONVOLUTION_METHODS[deconv_method]

coef_mat = deconv_function(overlap_matrix=overlap_matrix['mat'],am_areas = pixel_areas,cell_areas = cell_areas)

coef_mat_sparse = sparse.csr_matrix(coef_mat.values)

final_res = {'coef_mat' : coef_mat_sparse,
             'row_idx' : overlap_matrix['row_idx'],
             'col_idx' : overlap_matrix['col_idx']}

out_file_name = out_path + deconv_method + '.pkl'

with open(out_file_name, "wb") as f:
    pickle.dump(final_res, f)