require(sf)

source("../../../Utils/utils.R")

# Set to your own cluster storage path (matches Create_spatialdata_object.py's read of this file)
cluster_base_dir = Sys.getenv("SPATIAL_OMICS_CLUSTER_DIR", "/path/to/your/cluster/storage/Spatial_omics")

# Load data ---------------------------------------------------------------

atlas_geojson = sf::read_sf("../../../Data/Figure 1/Section_68_annots.geojson")
regions_res = readRDS("../../../Data/Figure 1/Final_regions_atlas.rds")
brain_ontology = read.csv("../../../Data/Figure 3/Brain_ontology_df.csv")

# Add missing from Cortex -----------------------------------------------------
region_annot = regions_res$regions_metadata
atlas_geojson$region[3] = "HPF"
atlas_geojson$id[3] = "1089"

missing_regions = setdiff(region_annot$acronym, atlas_geojson$region)

missing_gpd = atlas_geojson[1:length(missing_regions),]
for (i in 1:length(missing_regions)){
  reg_id = brain_ontology$id[brain_ontology$acronym == missing_regions[i]]
  reg_children = get_children_regions(brain_ontology, missing_regions[i])
  reg_children_geoms = atlas_geojson$geometry[atlas_geojson$region %in% reg_children]
  reg_polygon = sf::st_multipolygon(x = reg_children_geoms)

  missing_gpd$id[i] = reg_id
  missing_gpd$region[i] = missing_regions[i]
  missing_gpd$geometry[i] = reg_polygon
}

updated_atlas_geojson = rbind.data.frame(atlas_geojson, missing_gpd)

# Annotate fiber tracts ---------------------------------------------------
atlas_geojson = updated_atlas_geojson
affine_mat = rotation(180)
affine_mat[1,1] = 1

not_fiber_tract = atlas_geojson[atlas_geojson$region != "fiber tracts",]
fiber_tract_gpd = atlas_geojson[atlas_geojson$region == "fiber tracts",]

fiber_tract_gpd$region = paste0("fiber tracts_", c(1:nrow(fiber_tract_gpd)))

#Convert geometrycollection to polygon
for (i in 1:nrow(fiber_tract_gpd)){
  if (! "GEOMETRYCOLLECTION" %in% class(fiber_tract_gpd$geometry[[i]])){
    next()
  }
  else{
    fiber_tract_gpd$geometry[[i]] = fiber_tract_gpd$geometry[[i]] %>%
      st_union() %>%
      st_collection_extract("POLYGON") %>%
      st_cast("POLYGON")
  }
}

ft_parts = get_children_regions(brain_ontology, "fiber tracts")
ft_annot_parts = intersect(ft_parts, atlas_geojson$region)
ft_annot_geoms = atlas_geojson[atlas_geojson$region %in% ft_annot_parts,]

plot(fiber_tract_gpd$geometry * affine_mat, col = "darkgreen")
plot(st_combine(ft_annot_geoms$geometry) * affine_mat,
     add = T, col = '#ff333388')


#fiber tract sub-regions, identified by manual comparison against the Allen ontology reference

ft_annot_mapping = c(
  "fiber tracts_1" = "opt",
  "fiber tracts_2" = "mtt",
  "fiber tracts_3" = "fx",
  "fiber tracts_4" = "st",
  "fiber tracts_5" = "st",
  "fiber tracts_6" = "int",
  "fiber tracts_7" = "sm"
)

fiber_tract_gpd$region = str_replace_all(fiber_tract_gpd$region,
                                         ft_annot_mapping)
fiber_tract_gpd = fiber_tract_gpd[-c(8:10),]

final_atlas_gpd = rbind.data.frame(not_fiber_tract, fiber_tract_gpd)


# Save results ------------------------------------------------------------
saveRDS(final_atlas_gpd, "../../../Data/Figure 1/updated_atlas_geojson.rds")
sf::st_write(final_atlas_gpd, file.path(cluster_base_dir, "Data/SpatialData_input/updated_section_68_annots.geojson"),
             driver = "GeoJSON")
