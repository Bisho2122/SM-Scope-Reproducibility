library(tidyverse)
library(Seurat)
library(dbscan)
library(igraph)

# Data --------------------------------------------------------------------
metabo_transc = readRDS("../../../Data/Figure 1/Transc_metabo_seurat_combined_updated.rds")

# geopandas_intersect_cells_atlas.csv is produced by atlas_aggregate_job.py
# (Analysis_scripts/Upstream/Allen_Atlas_to_cells/Aggregate_atlas_cells_sdata.sh)
gpd_intersect_atlas = read.csv("../../../Data/Figure 1/geopandas_intersect_cells_atlas.csv")
gpd_intersect_atlas = gpd_intersect_atlas %>%
  dplyr::filter(region != "grey", cell_fraction > 0.95)

final_regions = readRDS("../../../Data/Figure 1/Final_regions_atlas.rds")

final_seurat_obj = metabo_transc

coord.df = final_seurat_obj@meta.data[,c("y","x")]

final_seurat_obj@images$image = new(Class = "SlideSeq",
                                    assay = "Spatial",
                                    key = "image_",
                                    coordinates = coord.df)


# Functions ---------------------------------------------------------------
convert_pixel_to_micrometers <- function(coords_df, pixel_size) {
  # Check if inputs are valid
  if (!is.data.frame(coords_df)) {
    stop("coords_df must be a data frame.")
  }
  if (!all(c("x", "y") %in% colnames(coords_df))) {
    stop("coords_df must contain columns named 'x' and 'y'.")
  }
  if (!is.numeric(pixel_size) || pixel_size <= 0) {
    stop("pixel_size must be a positive numeric value.")
  }

  # Convert coordinates
  coords_df$x <- coords_df$x * pixel_size
  coords_df$y <- coords_df$y * pixel_size

  # Return updated data frame
  return(coords_df)
}

calculate_snn_from_frnn <- function(frnn_result) {
  n <- length(frnn_result$id)  # Number of points
  snn_matrix <- matrix(0, nrow = n, ncol = n)  # Initialize SNN similarity matrix
  for (i in 1:n) {
    neighbors_i <- frnn_result$id[[i]]  # Neighbors of point i
    for (j in neighbors_i) {
      if (i != j) {
        neighbors_j <- frnn_result$id[[j]]  # Neighbors of point j
        # Calculate the number of shared neighbors
        shared_neighbors <- length(intersect(neighbors_i, neighbors_j))
        snn_matrix[i, j] <- shared_neighbors
      }
    }
  }

  return(snn_matrix)
}

Region_to_spatial_niches = function(reg_seurat_obj,
                                    radius_um = 100,
                                    min_cluster_size = 15,
                                    use_SNN_dist = T){
  reg_coords = reg_seurat_obj@meta.data[,c("x","y")] %>%
    convert_pixel_to_micrometers(pixel_size = 0.138)

  # snn_dist <- max(snn) - snn  # Transform to distance matrix

  if(!use_SNN_dist){
    dist = distances::distances(reg_coords) %>% as.matrix()
    dist[dist > radius_um] = 0
  }
  else{
    frnn = dbscan::frNN(reg_coords, eps = radius_um)

    snn = calculate_snn_from_frnn(frnn)
    rownames(snn) = rownames(reg_coords)
    colnames(snn) = rownames(reg_coords)

    dist = snn
  }

  graph <- igraph::graph_from_adjacency_matrix(dist,
                                       mode = "undirected",
                                       weighted = TRUE)

  louvain_clusters <- igraph::cluster_louvain(graph,resolution = 1)

  plot_data = reg_coords
  # Add cluster labels to the data
  plot_data$cluster <- factor(louvain_clusters$membership)

  sanity_plot = ggplot(plot_data, aes(x = x, y = y, color = cluster)) +
    geom_point(size = 2) +
    theme_minimal() +
    theme(legend.position="none")

  clust_sizes = igraph::sizes(louvain_clusters)
  clust_membs = igraph::membership(louvain_clusters)

  clust_to_consider = clust_sizes[clust_sizes >= min_cluster_size] %>%
    names()

  if(length(clust_to_consider) == 0){
    final_res = list("clusters" = NA,
                     "clust_res" = louvain_clusters,
                     "plot" = sanity_plot)
  }
  else{
    clust_membs = clust_membs[clust_membs %in% clust_to_consider]

    final_res = list("clusters" = data.frame(Cell = names(clust_membs),
                                             Niche = as.character(clust_membs)),
                     "clust_res" = louvain_clusters,
                     "plot" = sanity_plot)
  }

  return(final_res)


}

# Divide Regions into spatial niches - Spatial based only -----------------------------------------------------------------
n_cells_region = final_regions$final_gpd %>%
  dplyr::filter(X__index %in% colnames(final_seurat_obj)) %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(n_cells = n())

all_regions = unique(n_cells_region$region)


clust_results = list()
spatial_niches = list()
for (i in all_regions){

  reg_cells = final_regions$final_gpd %>%
    dplyr::filter(region == i) %>%
    dplyr::pull(X__index)
  reg_cells = intersect(reg_cells, colnames(final_seurat_obj))

  reg_seurat = SeuratObject:::subset.Seurat(final_seurat_obj,
                                            cells = reg_cells)

  clust_res = Region_to_spatial_niches(reg_seurat_obj = reg_seurat,
                                   radius_um = 100,
                                   min_cluster_size = 15,
                                   use_SNN_dist = T)
  clust_results[[i]] = clust_res
  if(!all(is.na(clust_res$clusters))){
    clust_res$clusters$Niche = paste0(i, "_", clust_res$clusters$Niche)

    clust_results[[i]] = clust_res
    spatial_niches[[i]] = clust_res$clusters

  }
}

spatial_niches = dplyr::bind_rows(spatial_niches, .id = "Region")

saveRDS(clust_results, "../../../Data/Figure 1/spatial_niches_all_clust_results.rds")
saveRDS(spatial_niches, "../../../Data/Figure 1/spatial_niches_regions.rds")

# Supp Tables 5 and 7 (cell-to-region assignment, spatial niches) are now generated by
# Analysis_scripts/Downstream/Generate_Supp_Tables.R from this script's saved Data/ objects.
