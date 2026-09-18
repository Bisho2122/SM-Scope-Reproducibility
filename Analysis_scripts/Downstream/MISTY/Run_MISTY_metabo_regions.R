library(tidyverse)
library(mistyR)
library(Seurat)


# Load data ---------------------------------------------------------------

metabo_RNA_seurat = readRDS("../../../Data/Figure 1/updated_Metabo_RNA_combined_seurat.rds")

final_regions = readRDS("../../../Data/Figure 1/Final_regions_atlas.rds")

# Prepare data for MISTY --------------------------------------------------
Region_composition = final_regions$final_gpd %>%
  dplyr::select(X__index, region) %>%
  dplyr::mutate(val = 1) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(names_from = region,
                     values_from = val,
                     values_fill = 0) %>%
  column_to_rownames("X__index")

colnames(Region_composition) = gsub("-", "_", colnames(Region_composition))

cells_metabo = as.matrix(GetAssayData(metabo_RNA_seurat,
                                      layer = "data",
                                      assay = "Metabo")) %>%
  t() %>%
  as.data.frame()
colnames(cells_metabo) = sub("[+-]", "_", colnames(cells_metabo))


# Run MISTY --------------------------------------------------

dir.create("../../../Data/Figure 3/updated_MISTY_metabo_regions")

out_dir = "../../../Data/Figure 3/updated_MISTY_metabo_regions"

common_cells = intersect(rownames(cells_metabo),
                         rownames(Region_composition))

comp_view = mistyR::create_initial_view(Region_composition[common_cells,])
colnames(comp_view$intraview$data) = gsub("-", "_", colnames(comp_view$intraview$data))

metabo_expr_view = mistyR::create_view("Metabo",
                               cells_metabo[common_cells,],
                               "Metabo")

combined_views = comp_view %>%
  mistyR::add_views(metabo_expr_view)

mistyR::run_misty(views = combined_views,results.folder = out_dir,bypass.intra = T)

misty_res = mistyR::collect_results(out_dir)
saveRDS(list("scenario_3" = list("Regions" = misty_res)), "../../../Data/Figure 3/updated_MISTY_metabo_region_result.rds")

misty_res %>%
  mistyR::plot_improvement_stats("intra.R2")%>%
  mistyR::plot_improvement_stats("gain.R2")

misty_res %>%
  mistyR::plot_view_contributions()

p = misty_res %>%
  mistyR::plot_interaction_heatmap("Metabo", clean = TRUE)
