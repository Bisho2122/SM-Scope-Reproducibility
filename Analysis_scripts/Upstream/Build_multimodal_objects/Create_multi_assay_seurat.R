library(Seurat)
library(tidyr)
library(tibble)
library(dplyr)

source("../../../Utils/utils.R")

# Load Data ---------------------------------------------------------------

transc_seurat = readRDS("../../../Data/Figure 1/Final_seurat_transc_updated.rds")
transc_seurat = transc_seurat$seurat_res$Seurat_obj

metabo_seurat = readRDS("../../../Data/Figure 1/updated_Metabo_Seurat_obj.rds")
colnames(metabo_seurat) = gsub("-", "_", colnames(metabo_seurat))

regions_res = readRDS("../../../Data/Figure 1/Final_regions_atlas.rds")

all_sg_res = readRDS("../../../Data/Figure 6/SG_eval_res_combined.rds")
spatial_niches = readRDS("../../../Data/Figure 1/spatial_niches_regions.rds")


# Add Cell type to transc metadata -----------------------------------------------
cell_type_meta = Idents(transc_seurat) %>%
  as.data.frame() %>%
  rownames_to_column("Cell")
colnames(cell_type_meta) = c("Cell","Cell_type_fine")

cell_type_meta$Cell_type_fine = str_replace_all(cell_type_meta$Cell_type_fine,
                                                c("Excit_CTX" = "Excitatory_VGLUT1",
                                                  "Inhib_CTX" = "Inhibitory_cortex",
                                                  "Excit_VGLUT2" = "Excitatory_VGLUT2",
                                                  "Excit_HPF" = "Excitatory_HPF",
                                                  "Inhib_Interneurons" = "Inhibitory_interneurons",
                                                  "Ventricular systems Glial cells" = "Ependymal cells"))

transc_seurat@meta.data = transc_seurat@meta.data %>%
  rownames_to_column("Cell") %>%
  dplyr::left_join(cell_type_meta) %>%
  column_to_rownames("Cell")

# Combining seurat objs ---------------------------------------------------
metabo_assay_name = names(metabo_seurat@assays)

metabo_RNA_seurat = transc_seurat

common_cells = Reduce(intersect,list(colnames(metabo_RNA_seurat),
                                     colnames(metabo_seurat)))

metabo_RNA_seurat = subset(metabo_RNA_seurat, cells = common_cells)

sc_metabo_M1_S1_S1 = subset_seurat_by_cells(metabo_seurat, common_cells)


metabo_RNA_seurat[["Metabo"]] = subset(sc_metabo_M1_S1_S1@assays$inverse_sampling_proportion_weighted_by_sampling_specificity_corrected,
                                cells = common_cells)

assertthat::are_equal(colnames(metabo_RNA_seurat@assays$RNA),
                      colnames(metabo_RNA_seurat@assays$Metabo))

metabo_RNA_seurat@images = sc_metabo_M1_S1_S1@images
metabo_RNA_seurat@images$image@coordinates$cells = gsub("-", "_",
                                                 metabo_RNA_seurat@images$image@coordinates$cells)

metabo_metadata = sc_metabo_M1_S1_S1@meta.data %>%
  dplyr::select(Cell,Area,seurat_clusters) %>%
  dplyr::distinct()
colnames(metabo_metadata) = c("Cell", "Cell_area", "Metabotype")
metabo_metadata$Metabotype = paste0("Met.", metabo_metadata$Metabotype)

all_metadata = metabo_RNA_seurat@meta.data %>%
  rownames_to_column("Cell") %>%
  dplyr::left_join(metabo_metadata) %>%
  column_to_rownames("Cell")

names(metabo_RNA_seurat@reductions) = paste0("RNA_", names(metabo_RNA_seurat@reductions))

for (i in c("pca", "umap")){
  metabo_RNA_seurat@reductions[[paste0("Metabo_", i)]] = sc_metabo_M1_S1_S1@reductions[[i]]
}

for (i in c("nn", "snn")){
  metabo_RNA_seurat@graphs[[paste0("Metabo_", i)]] = sc_metabo_M1_S1_S1@graphs[[paste0(metabo_assay_name, "_", i)]]
}


metabo_RNA_seurat@meta.data = all_metadata




# Add regions to metadata -------------------------------------------------
final_gpd = regions_res$final_gpd
region_annot = regions_res$regions_metadata
metadata = metabo_RNA_seurat@meta.data %>%
  rownames_to_column("Cell")
cell_to_region = metadata %>%
  dplyr::select(Cell) %>%
  dplyr::left_join(final_gpd[,c("X__index", "region")],
                   by = c("Cell" = "X__index")) %>%
  dplyr::filter(!is.na(region)) %>%
  dplyr::distinct() %>%
  dplyr::left_join(region_annot[,c("acronym", "name")],
                   by = c("region" = "acronym")) %>%
  dplyr::rename("region_name" = "name")

metabo_RNA_seurat@meta.data = metabo_RNA_seurat@meta.data %>%
  rownames_to_column("Cell") %>%
  dplyr::left_join(cell_to_region) %>%
  column_to_rownames("Cell")

metabo_RNA_seurat@meta.data$region[is.na(metabo_RNA_seurat@meta.data$region)] = "Uncertain"

# Add CT-Region to metadata -----------------------------------------------
metabo_RNA_seurat@meta.data$CT_Region = paste0(metabo_RNA_seurat@meta.data$Cell_type_fine,
                                               ".", metabo_RNA_seurat@meta.data$region)

# Add spatial niches to metadata ------------------------------------------
niches = spatial_niches %>% dplyr::select(Cell,Niche)
metabo_RNA_seurat@meta.data = metabo_RNA_seurat@meta.data %>%
  rownames_to_column("Cell") %>%
  dplyr::left_join(niches) %>%
  column_to_rownames("Cell")


# Add SpatialGlue Leiden to metadata  -------------------------------------------------
selected_param = "dim50_k15_metriceuclidean_wf5511.csv"

sel_sg_res = lapply(all_sg_res, function(x){
  x[x$param_comb == selected_param,]
})

sg_leiden = sel_sg_res$SG_clusts_metadata %>%
  dplyr::select(X, SG_leiden) %>%
  dplyr::distinct() %>%
  dplyr::mutate(SG_leiden = paste0("SG.", SG_leiden))

metabo_RNA_seurat@meta.data = metabo_RNA_seurat@meta.data %>%
  rownames_to_column("Cell") %>%
  dplyr::left_join(sg_leiden, by = c("Cell" = "X")) %>%
  column_to_rownames("Cell")


# Save seurat object ------------------------------------------------------

saveRDS(metabo_RNA_seurat, "../../../Data/Figure 1/updated_Metabo_RNA_combined_seurat.rds")
