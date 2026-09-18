import os
from metaspace import SMInstance
import metaspace
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from uuid import uuid4
import base64
import requests
from sm.engine.postprocessing.off_sample_wrapper import base64_images_to_doc, call_api
from cpyMSpec import isotopePattern, InstrumentModel

# Set to your own cluster storage path (matches Apply_QC_METASPACE_Annot.py's read/write of the same intermediates)
CLUSTER_BASE_DIR = os.environ.get("SPATIAL_OMICS_CLUSTER_DIR", "/path/to/your/cluster/storage/Spatial_omics")

def predict_off_sample(api_endpoint,aws_url):
    response = requests.get(aws_url)
    img_base64 = base64.b64encode(response.content).decode('utf-8')
    trial_img = base64_images_to_doc([img_base64])
    pred_doc = call_api(api_endpoint + '/predict', doc=trial_img)
    return pred_doc['predictions']
def calc_sparsity(arr_img):
    sparsity = np.count_nonzero(arr_img==0) / arr_img.size
    return sparsity

def Run_QC_pipeline(ds, db_name, db_version):
    ds_images = ds.all_annotation_images(only_first_isotope=True,fdr=1, database=(db_name, db_version))
    ann_docs = []
    for img in ds_images:
        formula = img._sf
        adduct = img._adduct
        ann_docs.append({
            'formula': formula,
            'adduct': adduct,
            'image': next(iter(img)),
            'image_url': img._urls[0],
        })
    off_sample_preds = [predict_off_sample(api_endpoint, d['image_url']) for d in tqdm(ann_docs, desc="Off sample pred")]
    all_sparsity_vals = [calc_sparsity(d['image']) for d in tqdm(ann_docs, desc="Calc Sparsity")]
    QC_df = []
    for i in range(len(off_sample_preds)):
        annot_res = {'formula' : ann_docs[i]['formula'],
                    'adduct' : ann_docs[i]['adduct'],
                    'off_sample_prob' : off_sample_preds[i][0]['prob'],
                    'off_sample_label' : off_sample_preds[i][0]['label'],
                    'Sparsity' : all_sparsity_vals[i]}
        QC_df.append(annot_res)
    QC_df = pd.DataFrame(QC_df)
    QC_df['DB_name'] = db_name
    QC_df['DB_version'] = db_version

    return (QC_df)


def _trim(mzs, ints, k):
    """Only keep top k peaks"""
    int_order = np.argsort(ints)[::-1]
    mzs = mzs[int_order][:k]
    ints = ints[int_order][:k]
    mz_order = np.argsort(mzs)
    mzs = mzs[mz_order]
    ints = ints[mz_order]
    return mzs, ints

def get_iso_theoritical(ds, ion_formula):

    iso_config = ds.config['isotope_generation']
    SIGMA_TO_FWHM = 2.3548200450309493  # 2 \sqrt{2 \log 2}

    iso_pattern = isotopePattern(str(ion_formula))
    iso_pattern.addCharge(int(iso_config['charge']))
    fwhm = float(iso_config['isocalc_sigma']) * SIGMA_TO_FWHM
    resolving_power = iso_pattern.masses[0] / fwhm
    instrument_model = InstrumentModel('tof', resolving_power)

    centr = iso_pattern.centroids(instrument_model)
    mzs_ = np.array(centr.masses)
    ints_ = 100.0 * np.array(centr.intensities)
    mzs_, ints_ = _trim(mzs_, ints_, iso_config['n_peaks'])

    n = len(mzs_)
    mzs = np.zeros(iso_config['n_peaks'])
    mzs[:n] = np.array(mzs_)
    ints = np.zeros(iso_config['n_peaks'])
    ints[:n] = ints_

    return mzs, ints

def calc_new_spectral_per_ion(iso_obs, iso_theo):
    assert isinstance(iso_obs, metaspace.sm_annotation_utils.IsotopeImages)

    iso_imgs_flat = np.array([img.flatten() for img in iso_obs._images])
    not_null = iso_imgs_flat[0] > 0
    obs_image_ints = [np.sum(iso_imgs_flat[i][not_null]) for i in range(iso_imgs_flat.shape[0])]

    theo_mz, theo_ints = iso_theo

    # assert np.all(np.around(theo_mz,2) == np.around(iso_obs._centroids, 2))

    obs_ints_norm = obs_image_ints / np.linalg.norm(obs_image_ints)
    theo_ints_norm = theo_ints / np.linalg.norm(theo_ints)
    theo_ints_norm[theo_ints_norm == 0] = 1e-10

    obs_theo_diff = abs(theo_ints_norm - obs_ints_norm)

    obs_theo_diff_norm = obs_theo_diff / np.maximum(obs_ints_norm,theo_ints_norm)
    obs_theo_diff_norm = 1 - obs_theo_diff_norm

    theo_weights_avg = theo_ints / sum(theo_ints)

    new_spectral = np.average(obs_theo_diff_norm, weights=theo_weights_avg)

    return new_spectral


sm = SMInstance()
api_endpoint = 'http://off-sample-api-load-balancer-630496755.eu-west-1.elb.amazonaws.com/off-sample'

ds_id = "2022-01-11_23h33m36s"
ds = sm.dataset(id = ds_id)

out_dir_path = f"{CLUSTER_BASE_DIR}/Results/Upstream/METASPACE_annotations/Annotation_QC"
all_dbs = [('HMDB', 'v4'), ('refmet', '27_06_2024'), ('CoreMetabolome', 'v3'),
           ('ChEBI', '2018-01'), ("RaMP_db_brain", "19_11_2024")]

for i in range(len(all_dbs)):
    database_name, database_version = all_dbs[i]
    db_QC = Run_QC_pipeline(ds, db_name = database_name, db_version = database_version)
    out_file_path = "/".join([out_dir_path, "_".join(['Annot_QC', ds_id ,database_name]) + '.csv'])
    db_QC.to_csv(out_file_path, index=False)
    
    
var_df = pd.read_csv(f"{CLUSTER_BASE_DIR}/Results/Upstream/Combined_obj_for_integration/M1_S1_S1_combined_adata_var_df.csv")

new_spectral_scores = []
for i in range(len(var_df.index)):
    iso_obs = ds.isotope_images(sf = var_df.formula[i], adduct=var_df.adduct[i])
    iso_theo = get_iso_theoritical(ds=ds, ion_formula=var_df.ionFormula[i])
    score = calc_new_spectral_per_ion(iso_obs=iso_obs, iso_theo=iso_theo)
    new_spectral_scores.append(score)


var_df['new_spectral'] = new_spectral_scores
var_df = var_df.set_index('ion_y')

Cl_annots = var_df[var_df.adduct == '+Cl']
Cl_annots['pass_spectral_Cl'] = True
Cl_annots['pass_spectral_Cl'][Cl_annots.new_spectral < 0.7] = False

var_df_with_cl_qc = pd.merge(left=var_df, right=Cl_annots['pass_spectral_Cl'], right_index=True, left_index=True, how='left')
var_df_with_cl_qc['pass_spectral_Cl'][var_df_with_cl_qc['pass_spectral_Cl'].isna()] = True

var_df_with_cl_qc.to_csv(f"{CLUSTER_BASE_DIR}/Results/Upstream/Combined_obj_for_integration/M1_S1_S1_combined_adata_var_df_w_new_spectral.csv",
              index=True)
              
