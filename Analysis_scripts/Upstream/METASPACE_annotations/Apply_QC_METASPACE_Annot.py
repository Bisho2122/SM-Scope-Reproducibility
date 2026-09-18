import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import tqdm as notebook_tqdm
#import squidpy as sq
import re

# Set to your own cluster storage path (matches METASPACE_Annot_QC.py's read/write of the same intermediates)
CLUSTER_BASE_DIR = os.environ.get("SPATIAL_OMICS_CLUSTER_DIR", "/path/to/your/cluster/storage/Spatial_omics")

def get_mz_range(mz,ppm):
    uncertainty_mz =  mz * (ppm * 1e-6)
    lower = mz - uncertainty_mz
    upper = mz + uncertainty_mz
    return lower, upper

def Map_isobars(var_df,ppm):

    # Cl_adduct_df = var_df[var_df['adduct'] == '+Cl']
    all_ions = var_df

    iso_bar_mapping = []

    for i in range(len(all_ions.index)):
        lower,upper = get_mz_range(all_ions['mz'][i], ppm=ppm)
        other_ions = var_df[var_df.mz.between(lower, upper)]
        other_ions = other_ions[other_ions['adduct'].isin(['-H'])]
        if len(other_ions.index) == 0:
            mapping = {'isobar' : all_ions.index[i],
                       'main_ion' : all_ions.index[i]}
        else:
            mapping = {'isobar' : list(other_ions.index),
                       'main_ion' : all_ions.index[i]}
        iso_bar_mapping.append(mapping)

    df = pd.DataFrame(iso_bar_mapping)

    if df['isobar'].apply(lambda x: isinstance(x, list)).any():
        df = df.explode('isobar', ignore_index=True)

    df_sorted = pd.DataFrame(
    {
        'isobar': df[['isobar', 'main_ion']].min(axis=1),
        'main_ion': df[['isobar', 'main_ion']].max(axis=1)
    })

    df = df_sorted
    df = df[df.main_ion != df.isobar]

    df = pd.merge(left=df, right=all_ions[['fdr','new_spectral']] ,left_on='isobar', right_index=True)
    df.rename(columns={'fdr': 'fdr_isobar',
                      'new_spectral' : 'new_spectral_isobar'}, inplace=True)
    df = pd.merge(left=df, right=all_ions[['fdr','new_spectral']] ,left_on='main_ion', right_index=True)
    df.rename(columns={'fdr': 'fdr_main_ion',
                      'new_spectral' : 'new_spectral_main_ion'}, inplace=True)

    return df

def filter_isobars(isobar_map_df):
    to_remove = []  # Initialize the result list

    for _, row in isobar_map_df.iterrows():

        if row['isobar'] in to_remove or row['main_ion'] in to_remove:
            continue

        # If C or D is NaN, compare E and F
        if pd.isna(row['fdr_isobar']) or pd.isna(row['fdr_main_ion']):
            if row['new_spectral_isobar'] < row['new_spectral_main_ion']:
                to_remove.append(row['isobar'])
            elif row['new_spectral_isobar'] > row['new_spectral_main_ion']:
                to_remove.append(row['main_ion'])
        # Otherwise, compare C and D
        elif row['fdr_main_ion'] < row['fdr_isobar']:
            to_remove.append(row['isobar'])
        elif row['fdr_main_ion'] > row['fdr_isobar']:
            to_remove.append(row['main_ion'])
        # If C == D, compare E and F
        else:
            if row['new_spectral_isobar'] < row['new_spectral_main_ion']:
                to_remove.append(row['isobar'])
            elif row['new_spectral_isobar'] > row['new_spectral_main_ion']:
                to_remove.append(row['main_ion'])

    return to_remove

def extract_elements(chemical_formula):
    # Use a regular expression to find elements (uppercase letter followed by optional lowercase letter)
    pattern = r"[A-Z][a-z]?"
    elements = re.findall(pattern, chemical_formula)
    # Return the unique elements
    return set(elements)


adata = sc.read_h5ad('../../../Data/Figure 1/M1_S1_S1_cell_ion_adata_w_multiple_deconv.h5ad')

adata.var.to_csv(f"{CLUSTER_BASE_DIR}/Results/Upstream/Combined_obj_for_integration/M1_S1_S1_combined_adata_var_df.csv")
updated_var_df = pd.read_csv(f"{CLUSTER_BASE_DIR}/Results/Upstream/Combined_obj_for_integration/M1_S1_S1_combined_adata_var_df_w_new_spectral.csv")
updated_var_df = updated_var_df.set_index('ion_y')

ions_unique = adata.var.drop_duplicates(subset=['ionFormula']).index
adata = adata[:,adata.var.index.isin(ions_unique)]

Isobar_ions = Map_isobars(updated_var_df, ppm=3)
to_remove_isobars = filter_isobars(Isobar_ions)

adata = adata[:,~adata.var.index.isin(to_remove_isobars)]

adata.var = pd.merge(left=adata.var, right=updated_var_df['pass_spectral_Cl'], left_index=True, right_index=True, how='left')
to_keep = adata.var[adata.var['pass_spectral_Cl']].index
adata = adata[:,adata.var.index.isin(to_keep)]

adata.write_h5ad('../../../Data/Figure 1/M1_S1_S1_cell_ion_adata_w_multiple_deconv.h5ad')
