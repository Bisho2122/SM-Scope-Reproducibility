import os
import spatialdata as sd
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, shape
from anndata import AnnData
import shapely

# Set to your own cluster storage path (matches Create_spatialdata_object.py's write of the same object)
CLUSTER_BASE_DIR = os.environ.get("SPATIAL_OMICS_CLUSTER_DIR", "/path/to/your/cluster/storage/Spatial_omics")

sdata = sd.SpatialData().read(f"{CLUSTER_BASE_DIR}/Results/Upstream/SpatialData/Final_SpatialData_obj_April_22_2025")

def aggregate(sdata_obj,ref_shape_name : str, ref_table_name : str, coord_system, ref_region_colname = "region",
              actual_val_key = "instance"):

    sdata_object = sdata_obj

    ref_table_shape = sdata_object.tables[ref_table_name].copy()
    old_region_name = ref_table_shape.obs.region_type.unique()[0]
    ref_table_shape.obs[ref_region_colname] = ref_table_shape.obs[ref_region_colname].str.replace(old_region_name,ref_shape_name)
    ref_table_shape.obs.instance = pd.to_numeric(ref_table_shape.obs.instance)
    del ref_table_shape.uns['spatialdata_attrs']
    sdata_object['table_shape'] = sd.models.TableModel.parse(ref_table_shape,region= ref_shape_name, region_key=ref_region_colname, instance_key='instance')

    values_ = sd._core.operations.aggregate._parse_element(element=ref_shape_name, sdata=sdata_object, str_for_exception="values")
    by_ = sd._core.operations.aggregate._parse_element(element="cell_polygons", sdata=sdata_object, str_for_exception="by")

    by_ = sd._core.operations.transform.transform(by_, to_coordinate_system= coord_system)
    values_ = sd._core.operations.transform.transform(values_, to_coordinate_system= coord_system)

    actual_values = sd._core.query.relational_query.get_values(
            value_key= actual_val_key, sdata=sdata_object, element_name= ref_shape_name, table_name= 'table_shape'
        )

    INDEX = "__index"

    value_key = [actual_val_key]
    for vk in value_key:
        if vk not in values_.columns:
            matched_values = actual_values[vk]
            matched_values.index = values_.index
            values_[vk] = matched_values

    by_[INDEX] = by_.index

    by_['cell_areas'] = by_.geometry.area

    overlayed = gpd.overlay(by_, values_, how='intersection', make_valid=False)

    df_1 = overlayed
    df_1['intersect_areas'] = overlayed.geometry.area
    df_1 = df_1[['__index', actual_val_key, 'intersect_areas']]
    df_2 = by_[['__index', 'cell_areas']]

    intersect_df = pd.merge(df_1, df_2, how='left')
    intersect_df['cell_fraction'] = intersect_df.intersect_areas / intersect_df.cell_areas

    return intersect_df


atlas_intersect_df = aggregate(sdata_obj = sdata, ref_shape_name = "atlas_annots", ref_table_name = "atlas_regions",
                               ref_region_colname = "region_type", coord_system = "atlas", actual_val_key = "region")
atlas_intersect_df.to_csv(f"{CLUSTER_BASE_DIR}/Results/Upstream/SpatialData/updated_geopandas_intersect_cells_atlas.csv")