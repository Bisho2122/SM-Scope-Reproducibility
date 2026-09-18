from allensdk.api.queries.ontologies_api import OntologiesApi
from allensdk.api.queries.image_download_api import ImageDownloadApi
from allensdk.api.queries.svg_api import SvgApi
from svgelements import *
import svgelements
from svgpath2mpl import parse_path
from shapely.geometry import Polygon
from shapely.validation import make_valid
import pandas as pd
import geopandas as gpd
import numpy as np

def flatten(nested_list):
    flat_list = []
    for item in nested_list:
        if isinstance(item, list):
            flat_list.extend(flatten(item))  # Recursively flatten the nested list
        else:
            flat_list.append(item)  # Append non-list items directly
    return flat_list

def brain_svg2geojson(svg_file,onto_df):
    parsed_svg = SVG.parse(svg_file)
    svg_elements = list(parsed_svg.elements())
    svg_elem_flat = flatten(svg_elements)
    brain_polygons = []
    brain_ids = []
    region_names = []
    for p in range(len(svg_elem_flat)):
        svg_path = svg_elem_flat[p]
        onto_id = svg_path.values['structure_id']
        onto_name = onto_df.acronym[onto_df.id == int(onto_id)].values[0]

        svg_pts = [i for i in svg_path.as_points()]
        poly_intermediate = svgelements.Polygon(svg_pts)
        mpl_path = parse_path(poly_intermediate.d())
        coordinates = mpl_path.to_polygons()
        poly = Polygon(coordinates[0])

        if not poly.is_valid:
            poly = make_valid(poly)

        brain_polygons.append(poly)
        brain_ids.append(onto_id)
        region_names.append(onto_name)
    df = {'id': brain_ids, 'region': region_names, 'geometry': brain_polygons}
    gdf = gpd.GeoDataFrame(df)
    gdf = gdf.drop_duplicates(subset='geometry')
    return gdf


atlas_id = 1
atlas_struct_id = 1
atlas_struct_name = 'Mouse Brain Atlas'
data_export_path  = "../../Data/Affinder_atlas/Allen_atlas_sections/"

image_api = ImageDownloadApi()
image_sections = image_api.atlas_image_query(atlas_id=atlas_id)

onto = OntologiesApi()
main_onto = onto.get_structures_with_sets(structure_graph_ids=atlas_struct_id,count=True)

onto_df = pd.DataFrame(main_onto)
onto_df.to_csv('../../../Data/Figure 3/Brain_ontology_df.csv')

for i in range(66, 75):
    sec_name = "Section_" + str(i)
    sec_svg_out = data_export_path + sec_name + '.svg'

    image_id = image_sections[i]['id']

    for annot in [True, False]:
        if annot :
            sec_jpg_out = data_export_path + sec_name + '_annotated.jpg'
        else:
            sec_jpg_out = data_export_path + sec_name + '.jpg'
        image_api.download_atlas_image(atlas_image_id=image_id,annotation = annot, projection = False, atlas = 1,
                         downsample = 0, file_path = sec_jpg_out)

    SvgApi().download_svg(section_image_id=image_id, file_path=sec_svg_out)


svg_input_path = "../../Data/Affinder_atlas/Allen_atlas_sections/Section_68.svg"
gdf_68 = brain_svg2geojson(svg_input_path,onto_df)
gdf_68.to_file("../../Results/Affinder_atlas/Section_68_annots.geojson", driver='GeoJSON')