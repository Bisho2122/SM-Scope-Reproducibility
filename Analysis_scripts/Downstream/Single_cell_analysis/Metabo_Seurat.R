require(Seurat)
require(anndataR)
library(tidyverse)

# Load Data ---------------------------------------------------------------
multi_layer_deconv_adata = anndataR::read_h5ad("../../../Data/Figure 1/M1_S1_S1_cell_ion_adata_w_multiple_deconv.h5ad")

multi_layer_deconv_seurat = anndataR::to_Seurat(multi_layer_deconv_adata)

manual_annot_flags = read.csv("../../../Data/General/updated_annots_manual_curation.csv",
                              na.strings = c("", "NA"))


# Feature selection -------------------------------------------------------
feats_to_keep = manual_annot_flags %>%
  dplyr::filter(FP_flag != 1)
feats_to_keep$ion[!is.na(feats_to_keep$possible_isobar)] = feats_to_keep$possible_isobar[!is.na(feats_to_keep$possible_isobar)]
feats_to_keep = feats_to_keep$ion[!duplicated(feats_to_keep$ion)]

# Prepare Seurat object ----------------------------------------------------------
assay_name = "inverse_sampling_proportion_weighted_by_sampling_specificity_corrected"
counts_mat = multi_layer_deconv_adata$layers[[assay_name]] %>% t()
rownames(counts_mat) = multi_layer_deconv_adata$var_names
colnames(counts_mat) = multi_layer_deconv_adata$obs_names

metabo_seurat = CreateSeuratObject(
  counts = counts_mat,
  assay = assay_name,
  meta.data = multi_layer_deconv_seurat@meta.data
)

metabo_seurat = metabo_seurat[feats_to_keep,]

coord.df = metabo_seurat@meta.data[,c("y","x")]
coord.df$y = -1 * coord.df$y
coord.df$x = -1 * coord.df$x


metabo_seurat@images$image = new(Class = "SlideSeq",
                                 assay = "Spatial",
                                 key = "image_",
                                 coordinates = coord.df)


# Run Seurat --------------------------------------------------------------
metabo_seurat = NormalizeData(metabo_seurat,
                              normalization.method = "LogNormalize",
                              scale.factor = 10000)

metabo_seurat <- ScaleData(metabo_seurat, do.scale = T,
                           do.center = T)
var_feats = rownames(metabo_seurat)
metabo_seurat <- RunPCA(object = metabo_seurat, approx = F,
                        features = var_feats,
                        assay = assay_name)

metabo_seurat <- FindNeighbors(metabo_seurat,
                               dims = 1:15, k.param = 50)

metabo_seurat <- FindClusters(metabo_seurat,
                              resolution = 0.8,
                              algorithm = 4)

# clustree_data = clustree(metabo_seurat, prefix = paste0(assay_name, "_snn_res."),
#                          return = "graph") %>% igraph::as_data_frame()

metabo_seurat <- RunUMAP(metabo_seurat,
                         dims = 1:15,
                         reduction = "pca",
                         reduction.name = "umap",
                         min.dist = 0.01, spread = 2, metric = "euclidean")

saveRDS(metabo_seurat, "../../../Data/Figure 1/updated_Metabo_Seurat_obj.rds")
