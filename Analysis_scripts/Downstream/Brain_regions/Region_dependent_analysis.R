library(tidyverse)

# load_data ---------------------------------------------------------------

# geopandas_intersect_cells_atlas.csv is produced by atlas_aggregate_job.py
# (Analysis_scripts/Upstream/Allen_Atlas_to_cells/Aggregate_atlas_cells_sdata.sh)
gpd_intersect_atlas = read.csv("../../../Data/Figure 1/geopandas_intersect_cells_atlas.csv")
gpd_intersect_atlas = gpd_intersect_atlas %>%
  dplyr::filter(region != "grey", cell_fraction > 0.95)

brain_ontology = read.csv("../../../Data/Figure 3/Brain_ontology_df.csv")

# Function ----------------------------------------------------------------

get_children_regions = function(ontology_df, ROI_acronym){
  id = ontology_df$id[ontology_df$acronym == ROI_acronym]
  child_paths = ontology_df$structure_id_path[str_detect(ontology_df$structure_id_path,
                                                         paste0(id, "/"))]
  children = sub(paste0(".*", id), "", child_paths)
  children = strsplit(children, "/") %>%
    unlist() %>%
    unique()
  children = children[children != ""] %>% as.integer()
  children_regions = ontology_df$acronym[ontology_df$id %in% children]
  return(children_regions)
}
make_exclusive_regions = function(gpd, ontology_df,
                                  min_cells_region = 15,
                                  keep_regions = NULL){
  updated_gpd = gpd
  cells_per_region = gpd %>%
    dplyr::group_by(region) %>%
    dplyr::summarise(n_cells = n())

  final_regions = ontology_df %>%
    dplyr::filter(acronym %in% gpd$region) %>%
    dplyr::select(acronym, depth, name, id, parent_structure_id,
                  structure_id_path) %>%
    dplyr::left_join(cells_per_region, by = c("acronym" = "region"))

  if(!is.null(keep_regions)){
    to_keep = keep_regions
  }
  else{
    to_keep = c()
  }

  all_regions = unique(final_gpd$region)
  for (r in all_regions){
    r_children = get_children_regions(ontology_df = final_regions,
                                      ROI_acronym = r)
    if(length(r_children) == 0 || r %in% to_keep){
      next()
    }
    else{
      remove_children = final_regions %>%
        dplyr::filter(acronym %in% r_children, n_cells < min_cells_region) %>%
        pull(acronym)
      if(length(remove_children) == length(r_children)){
        to_remove = remove_children
      }
      else{
        to_remove = c(r,remove_children)
      }

      #Updating gpd
      updated_gpd = updated_gpd[!updated_gpd$region %in% to_remove,]
    }
  }
  return(updated_gpd)
}
remove_cells_mixed_regions = function(gpd, min_cells_region = 15){
  #Check mutual exclusiveness
  regions_per_cell = gpd %>%
    dplyr::group_by(X__index) %>%
    dplyr::summarise(n_regions = length(unique(region)))

  cells_mixed_regions = regions_per_cell$X__index[regions_per_cell$n_regions != 1]

  updated_gpd = gpd[!gpd$X__index %in% cells_mixed_regions,]

  cells_per_region = updated_gpd %>%
    dplyr::group_by(region) %>%
    dplyr::summarise(n_cells = n())

  to_remove_region = cells_per_region$region[cells_per_region$n_cells < min_cells_region]

  updated_gpd = updated_gpd[!updated_gpd$region %in% to_remove_region,]

  return(updated_gpd)
}

# Prepare final non-overlapping regions ---------------------------------------------------

#Regions were visually inspected against online reference and the
#following regions are added that should be mutually exclusive by shape


final_gpd = gpd_intersect_atlas
no_depth_beyond = c("Isocortex","OLF", "RHP", "CTXsp")
include_regions = c("cc","ZI")

for (r in no_depth_beyond){
  x = get_children_regions(ontology_df = brain_ontology, ROI_acronym = r)
  final_gpd$region[final_gpd$region %in% x] = r
}

final_gpd = make_exclusive_regions(gpd = final_gpd,
                                   ontology_df = brain_ontology,
                                   keep_regions = include_regions)
final_gpd = remove_cells_mixed_regions(final_gpd)

cells_per_region = final_gpd %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(n_cells = n())

final_regions = brain_ontology %>%
  dplyr::filter(acronym %in% final_gpd$region) %>%
  dplyr::select(acronym, depth, name, id, parent_structure_id,
                structure_id_path) %>%
  dplyr::left_join(cells_per_region, by = c("acronym" = "region"))

#Check mutual exclusiveness
regions_per_cell = final_gpd %>%
  dplyr::group_by(X__index) %>%
  dplyr::summarise(n_regions = length(unique(region)))

saveRDS(list("final_gpd" = final_gpd,
             "regions_metadata" = final_regions),
        "../../../Data/Figure 1/Final_regions_atlas.rds")
