import os
import spatialdata as sd
from imageio.v3 import imread
import anndata as ad
from metaspace_converter import metaspace_to_anndata
import numpy as np
from metaspace import SMInstance
import pandas as pd
from skimage import measure
import geopandas as gpd
from shapely.geometry import Polygon, shape, MultiPolygon, GeometryCollection

from PIL import Image

# Set to your own cluster storage path
CLUSTER_BASE_DIR = os.environ.get("SPATIAL_OMICS_CLUSTER_DIR", "/path/to/your/cluster/storage/Spatial_omics")


"""
Functions
"""

def create_square_polygon_vectorized(df, side_length=0.5):
    """
    Create a square polygon from the given X and Y coordinates with a specified side length.

    Parameters:
        x (float or numpy.ndarray): X coordinate(s) of the square's bottom-left corner.
        y (float or numpy.ndarray): Y coordinate(s) of the square's bottom-left corner.
        side_length (float): Length of each side of the square.

    Returns:
        geopandas.GeoDataFrame: GeoDataFrame containing the square polygon.
    """
    # Convert inputs to numpy arrays if not already
    x_centroids = df['X']
    y_centroids = df['Y']

    bottom_left_x = x_centroids - side_length
    bottom_left_y = y_centroids - side_length
    bottom_right_x = x_centroids + side_length
    bottom_right_y = y_centroids - side_length
    top_right_x = x_centroids + side_length
    top_right_y = y_centroids + side_length
    top_left_x = x_centroids - side_length
    top_left_y = y_centroids + side_length

    # Create a DataFrame with the calculated coordinates
    polygons_df = pd.DataFrame({
        'bottom_left_x': bottom_left_x,
        'bottom_left_y': bottom_left_y,
        'bottom_right_x': bottom_right_x,
        'bottom_right_y': bottom_right_y,
        'top_right_x': top_right_x,
        'top_right_y': top_right_y,
        'top_left_x': top_left_x,
        'top_left_y': top_left_y
    })

    # Create Polygon objects
    polygons = [
        Polygon(((row['bottom_left_x'], row['bottom_left_y']),
                    (row['bottom_right_x'], row['bottom_right_y']),
                    (row['top_right_x'], row['top_right_y']),
                    (row['top_left_x'], row['top_left_y']),
                    (row['bottom_left_x'], row['bottom_left_y'])))
        for index, row in polygons_df.iterrows()
    ]

    # Convert the list of Polygons to a GeoDataFrame
    gdf = gpd.GeoDataFrame(geometry=polygons)

    return gdf

def flatten_geometry_collection(geom):
    if isinstance(geom, GeometryCollection):
        geometries = [part for part in geom.geoms if part.geom_type in ['Polygon', 'MultiPolygon']]
        if len(geometries) == 1:
            return geometries[0]  # Return single Polygon or MultiPolygon directly
        elif len(geometries) > 1:
            return MultiPolygon(geometries)  # Return MultiPolygon if multiple geometries
        else:
            return None  # Return None if no valid geometries found
    elif geom.geom_type in ['Polygon', 'MultiPolygon']:
        return geom
    else:
        return None  # Handle other geometry types as needed


"""
Read and import elements
"""

#Read images
dapi_arr = imread(f"{CLUSTER_BASE_DIR}/Data/SpatialData_input/32753-Slide1_A1_DAPI.tiff")
dapi_arr = np.expand_dims(dapi_arr,axis=0) #To add a dimension for the channel

#Read detections
cell_arr = imread(f"{CLUSTER_BASE_DIR}/Results/QuPath/32753-Slide1_A1_DAPI-detections.tif")

#Read BF image
bf_arr = imread(f"{CLUSTER_BASE_DIR}/Data/SpatialData_input/M1_Sand1_Sec1_cropped.tif")
bf_arr = np.expand_dims(bf_arr,axis=0) #To add a dimension for the channel

#Read ion image tif
ion_image_arr = imread(f"{CLUSTER_BASE_DIR}/Data/SpatialData_input/M1_Sand1_sec1_ion.tif")
#ion_image_arr = np.expand_dims(ion_image_arr,axis=0) #To add a dimension for the channel

# Read annotated atlas image
atlas_img = imread(f"{CLUSTER_BASE_DIR}/Data/SpatialData_input/Section_68_annotated.jpg")

#Read transformations
BF_cell = np.loadtxt(f"{CLUSTER_BASE_DIR}/Data/SpatialData_input/M1_Sand1_Sec1_to_ISS.txt",delimiter=',')
ion_BF = np.loadtxt(f"{CLUSTER_BASE_DIR}/Data/SpatialData_input/Ion_to_M1_Sand1_Sec1.txt",delimiter=',')
atlas_cell = np.loadtxt(f"{CLUSTER_BASE_DIR}/Data/SpatialData_input/Final_section_68_whole_to_cells_HPF_mat.txt",delimiter=',')

#Read Resolve adata
adata = ad.read_h5ad('../../../Data/Figure 1/Resolve_adata.h5ad')
adata.obs['region'] = 'cells'
adata.obs['instance'] = adata.obs.index
adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})

updated_instances_values = np.array([int(i.split("Cell_")[1]) for i in adata.obs.instance.values])
adata.obs.instance = pd.Series(updated_instances_values, index=adata.obs.index)

"""
Create Metabo adata
"""
# Create Metabo adata
sm = SMInstance()
M1_S1_S1_ds_id= '2022-01-11_23h33m36s'
metabo_adata = metaspace_to_anndata(dataset_id=M1_S1_S1_ds_id, fdr = 0.1,database=('CoreMetabolome','v3'),
                                    sm=sm,)

flattened_ion_labels = metabo_adata.obs.index.to_numpy()
metabo_adata.obs['region'] = 'ions'
metabo_adata.obs['instance'] = flattened_ion_labels.astype(np.int32)

"""
Create Polygons
"""

#Create polygons from cell detections
indiv_cells = measure.regionprops(cell_arr.T)
cell_num = pd.Series({"Cell_" + str(indiv_cells[k].label) : k for k in range(len(indiv_cells))})
gdf_cells = gpd.GeoDataFrame(
        index=cell_num.index, geometry=[Polygon(indiv_cells[cell_num].coords) for cell_num in cell_num]
    )
gdf_cells = gdf_cells[gdf_cells.index.isin(adata.obs.index)].reindex(adata.obs.index)

#Ion image polygons
ion_spatial_coords = pd.DataFrame(metabo_adata.obsm['spatial'])
ion_spatial_coords.columns = ['X', 'Y']
ion_spatial_coords['X'] = ion_spatial_coords['X'] + 0.5
ion_spatial_coords['Y'] = ion_spatial_coords['Y'] + 0.5
ion_image_polygon = create_square_polygon_vectorized(ion_spatial_coords, side_length=0.5)

# Atlas polygons

to_remove = ["Isocortex","OLF", "RHP", "CTXsp"]

atlas_gdf = gpd.read_file(f"{CLUSTER_BASE_DIR}/Data/SpatialData_input/updated_section_68_annots.geojson")
atlas_gdf['geometry'] = atlas_gdf['geometry'].apply(flatten_geometry_collection)
atlas_gdf = atlas_gdf[~atlas_gdf.region.isin(to_remove)]

atlas_gdf_table = pd.DataFrame(atlas_gdf.drop(columns='geometry'))
atlas_gdf = atlas_gdf.drop(columns=['id','region'])

"""
Parsing spatialdata elements
"""

cell_polygons = sd.models.ShapesModel.parse(gdf_cells,
	transformations = {'RNA' : sd.transformations.Identity()})

dapi_img = sd.models.Image2DModel.parse(dapi_arr, dims = ('c', 'y','x'), c_coords = ['dapi'],
	transformations = {'RNA' : sd.transformations.Identity()})

cell_labels = sd.models.Labels2DModel.parse(cell_arr, dims = ('y','x'),
	transformations = {'RNA' : sd.transformations.Identity()})

bf_img = sd.models.Image2DModel.parse(bf_arr, dims = ('c', 'y','x'), c_coords = 'BF',
	transformations = {'RNA' : sd.transformations.Affine(BF_cell,
		input_axes = ('y','x'),
		output_axes = ('y','x'))})

ion_x = metabo_adata.obs.ion_image_shape_x.max()
ion_y = metabo_adata.obs.ion_image_shape_y.max()

ion_labels_arr = flattened_ion_labels.reshape((ion_y,ion_x))
ion_labels_arr = ion_labels_arr.astype(np.uint32)
#ion_labels_arr = np.transpose(ion_labels_arr)

ion_labels = sd.models.Labels2DModel.parse(ion_labels_arr, dims = ('y','x'),
	transformations = {'Metabo' : sd.transformations.Affine(ion_BF,
		input_axes = ('y','x'),
		output_axes = ('y','x'))})

ion_image = sd.models.Image2DModel.parse(ion_image_arr, dims = ('y','x','c'), c_coords = ['R','G','B'],
	transformations = {'Metabo' : sd.transformations.Affine(ion_BF,
		input_axes = ('y','x'),
		output_axes = ('y','x'))}) #Check input axes and output axes

ion_image_shape = sd.models.ShapesModel.parse(ion_image_polygon,
	transformations = {'Metabo' : sd.transformations.Affine(ion_BF,
		input_axes = ('y','x'),
		output_axes = ('y','x'))})

atlas_image = sd.models.Image2DModel.parse(atlas_img, dims = ('y','x','c'), c_coords = ['R','G','B'],
	transformations = {'RNA' : sd.transformations.Affine(atlas_cell,
		input_axes = ('y','x'),
		output_axes = ('y','x'))}) #Check input axes and output axes

atlas_polygons = sd.models.ShapesModel.parse(atlas_gdf,
	transformations = {'RNA' : sd.transformations.Affine(atlas_cell,
		input_axes = ('y','x'),
		output_axes = ('y','x'))})

"""
Create Multi-table adata
"""
cells_table = sd.models.TableModel.parse(adata, region_key = 'region', instance_key = 'instance',
	region = 'cells')

ions_table = sd.models.TableModel.parse(metabo_adata, region_key = 'region', instance_key = 'instance',
	region = 'ions')

atlas_adata = ad.AnnData(obs = atlas_gdf_table)
atlas_adata.obs['region_type'] = "atlas"
atlas_adata.obs['instance'] = atlas_adata.obs.index
atlas_table = sd.models.TableModel.parse(atlas_adata, region_key = 'region_type', instance_key = 'instance',
	region = 'atlas')

"""
Create SpatialData object
"""

#Create SpatialData object
sdata = sd.SpatialData(
	images = {'dapi' : dapi_img,
	'bf_maldi' : bf_img,
    'ion_image' : ion_image,
    'atlas_section' : atlas_image},
	labels = {'cells' : cell_labels,
	'ions' : ion_labels},
    shapes = {'cell_polygons' : cell_polygons,
              'ion_image_polyg': ion_image_shape,
              'atlas_annots' : atlas_polygons},
	table = {'cell_table' : cells_table,
          'ion_table' : ions_table})
sdata.tables['atlas_regions'] = atlas_table
sdata.shapes['atlas_annots'] = atlas_polygons
sdata.tables['atlas_regions'] = atlas_table

"""
Apply transformations and align coordinate systems
"""
## Get transformations
moving_bf = sdata.images['bf_maldi']
moving_ions = sdata.labels['ions']
moving_ion_image = sdata.images['ion_image']
moving_ion_polyg = sdata.shapes['ion_image_polyg']
moving_atlas_image = sdata.images['atlas_section']
moving_atlas_polygons = sdata.shapes['atlas_annots']

moving_coord_system = 'Metabo'

ref_RNA = sdata.labels['cells']
ref_RNA_polyg = sdata.shapes['cell_polygons']
ref_coord_system = 'RNA'

#sd.transformations.Sequence
affine = sd.transformations.get_transformation(moving_bf, ref_coord_system)

atlas_affine = sd.transformations.get_transformation(moving_atlas_image, ref_coord_system)

metabo_moving_transformation = sd.transformations.get_transformation(moving_ions,
                                                                     moving_coord_system)
RNA_ref_transformation = sd.transformations.get_transformation(ref_RNA,
                                                                     ref_coord_system)
assert isinstance(metabo_moving_transformation, sd.transformations.BaseTransformation)
assert isinstance(RNA_ref_transformation, sd.transformations.BaseTransformation)

#Apply transformation
new_coord_sys = "aligned"

new_moving_transformation = sd.transformations.Sequence([metabo_moving_transformation, affine])
new_reference_transformation = RNA_ref_transformation

sd.transformations.set_transformation(moving_ions, new_moving_transformation, new_coord_sys)
sd.transformations.set_transformation(moving_ion_image, new_moving_transformation, new_coord_sys)
sd.transformations.set_transformation(moving_ion_polyg, new_moving_transformation, new_coord_sys)

sd.transformations.set_transformation(ref_RNA, new_reference_transformation, new_coord_sys)
sd.transformations.set_transformation(ref_RNA_polyg, new_reference_transformation, new_coord_sys)

#Apply transformation
new_coord_sys = "atlas"

sd.transformations.set_transformation(moving_atlas_image, atlas_affine, new_coord_sys)
sd.transformations.set_transformation(moving_atlas_polygons, atlas_affine, new_coord_sys)

sd.transformations.set_transformation(ref_RNA, new_reference_transformation, new_coord_sys)
sd.transformations.set_transformation(ref_RNA_polyg, new_reference_transformation, new_coord_sys)

"""
Save Spatialdata object
"""

sdata.write(f"{CLUSTER_BASE_DIR}/Results/Upstream/SpatialData/Final_SpatialData_obj_April_22_2025",overwrite=True)