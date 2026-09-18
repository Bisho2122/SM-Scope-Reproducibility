library(anndataR)
library(Seurat)
library(tidyverse)

transc_seurat = readRDS("../../../Data/Figure 1/Final_seurat_transc_updated.rds")

counts = transc_seurat$seurat_res$Seurat_obj@assays$RNA@counts
spatial_coords = transc_seurat$spatial_coords %>%
  dplyr::select(Cell, Area, x, y) %>%
  dplyr::distinct()
metadata = transc_seurat$seurat_res$Seurat_obj@meta.data %>%
  rownames_to_column("Cell") %>%
  dplyr::left_join(spatial_coords) %>%
  column_to_rownames("Cell")


all_genes = rownames(counts)

gene_metadata = data.frame(feature = all_genes)
rownames(gene_metadata) = all_genes

obj = SeuratObject::CreateSeuratObject(counts = counts,meta.data = metadata)
obj[["RNA"]] = AddMetaData(object = obj[["RNA"]], gene_metadata)

adata = as_AnnData(obj,x_mapping = "counts")
adata$write_h5ad(path = "../../../Data/Figure 1/Resolve_adata.h5ad",mode = "w")
