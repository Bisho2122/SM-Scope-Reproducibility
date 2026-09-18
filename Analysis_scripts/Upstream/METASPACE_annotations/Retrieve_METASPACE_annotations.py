import os
from metaspace import SMInstance
import numpy as np
from pandas import json_normalize
import pandas as pd

# Set to your own cluster storage path
CLUSTER_BASE_DIR = os.environ.get("SPATIAL_OMICS_CLUSTER_DIR", "/path/to/your/cluster/storage/Spatial_omics")

sm = SMInstance()
ds = sm.dataset(id='2022-01-11_23h33m36s')

ds_id = '2022-01-11_23h33m36s'
out_dir = f"{CLUSTER_BASE_DIR}/Results/Upstream/METASPACE_annotations/"

db_res = []
for db in ds.database_details:
    annot_filter = {'fdrLevel' : 1.0,
            'hasChemMod': False,
            'hasNeutralLoss': False,
            'databaseId' : db.id
            }
    ds = sm.dataset(id = ds_id)
    records = ds._gqclient.getAnnotations(
            annotationFilter=annot_filter,
            datasetFilter={'ids': ds.id},
            colocFilter=None,
        )
    df = json_normalize(records)
    df['database'] = db.name
    db_res.append(df)

final = pd.concat(db_res)
out_path = out_dir + ds_id + '.csv'
final.to_csv(out_path)
