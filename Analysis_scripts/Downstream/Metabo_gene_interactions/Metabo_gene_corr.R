library(tidyverse)
library(Seurat)
library(scLink)
library(ComplexHeatmap)

source("../../../Utils/utils.R")


# Load Data ---------------------------------------------------------------

metabo_RNA_seurat = readRDS("../../../Data/Figure 1/updated_Metabo_RNA_combined_seurat.rds")


# Correlation spatial niches ----------------------------------------------
all_niches = metabo_RNA_seurat@meta.data$Niche %>%
  unique() %>% na.omit()

pb = txtProgressBar(min = 0, max = length(all_niches), style = 3)

all_niches_corr = list()
for (i in 1:length(all_niches)){
  all_niches_corr[[all_niches[i]]] = calc_rna_metabo_corr(seurat_obj = metabo_RNA_seurat,
                                              metadata_filters = list("Niche" = all_niches[i]),
                                              corr_method = "spearman",
                                              min_nonzero_cells = 30,
                                              feat_nonzero_prop = 0.5,
                                              metabo_LOD = 0,
                                              nonzero_only = F)
  setTxtProgressBar(pb, i)
}

niches_corr = lapply(all_niches_corr, function(x){x$cor_df}) %>%
  dplyr::bind_rows(.id = "Niche") %>%
  dplyr::filter(feat_A_type != feat_B_type,
                feat_A_type == "RNA",
                feat_B_type == "Metabo")

niches_corr_input = lapply(all_niches_corr, function(x){x$cor_input})
names(niches_corr_input) = names(all_niches_corr)

saveRDS(niches_corr,"../../../Data/Figure 4/Per_niche_gene_metabo_corr.rds")


# Correlation SG leiden ---------------------------------------------------
all_sg = metabo_RNA_seurat@meta.data$SG_leiden %>%
  unique() %>% na.omit()

pb = txtProgressBar(min = 0, max = length(all_sg), style = 3)

all_sg_corr = list()
for (i in 1:length(all_sg)){
  all_sg_corr[[all_sg[i]]] = calc_rna_metabo_corr(seurat_obj = metabo_RNA_seurat,
                                                          metadata_filters = list("SG_leiden" = all_sg[i]),
                                                          corr_method = "spearman",
                                                          min_nonzero_cells = 30,
                                                          feat_nonzero_prop = 0.5,
                                                          metabo_LOD = 0,
                                                          nonzero_only = F)
  setTxtProgressBar(pb, i)
}

sg_corr_df = lapply(all_sg_corr, function(x){x$cor_df}) %>%
  dplyr::bind_rows(.id = "SG_leiden") %>%
  dplyr::filter(feat_A_type != feat_B_type,
                feat_A_type == "RNA",
                feat_B_type == "Metabo")

sg_corr_df = sg_corr_df %>% dplyr::filter(n_cells_cor >= 30)

sg_corr_input = lapply(all_sg_corr, function(x){x$cor_input})
names(sg_corr_input) = names(all_sg_corr)

sg_corr_top_per_pair = sg_corr_df %>%
  dplyr::group_by(feat_A, feat_B) %>%
  dplyr::slice_max(abs(coef), n = 1, with_ties = F)

length(which(abs(sg_corr_top_per_pair$coef) >= 0.5))


# Correlation all cells ---------------------------------------------------

all_cells_corr = calc_rna_metabo_corr(seurat_obj = metabo_RNA_seurat,
                                      metadata_filters = list(),
                                      corr_method = "spearman",
                                      min_nonzero_cells = 30,
                                      feat_nonzero_prop = 0.9,
                                      metabo_LOD = 0,
                                      nonzero_only = F)

all_cells_corr_df = all_cells_corr$cor_df %>%
  dplyr::filter(feat_A_type != feat_B_type,
                feat_A_type == "RNA",
                feat_B_type == "Metabo")

all_cells_corr_input = all_cells_corr$cor_input

# Correlation per region --------------------------------------------------

all_regions = metabo_RNA_seurat@meta.data$region %>%
  unique() %>% na.omit()
all_regions = all_regions[all_regions != "Uncertain"]

pb = txtProgressBar(min = 0, max = length(all_regions), style = 3)

all_reg_corr = list()
for (i in 1:length(all_regions)){
  all_reg_corr[[all_regions[i]]] = calc_rna_metabo_corr(seurat_obj = metabo_RNA_seurat,
                                                          metadata_filters = list("region" = all_regions[i]),
                                                          corr_method = "spearman",
                                                          min_nonzero_cells = 30,
                                                          feat_nonzero_prop = 0.5,
                                                          metabo_LOD = 0,
                                                          nonzero_only = F)
  setTxtProgressBar(pb, i)
}

all_reg_corr_df = lapply(all_reg_corr, function(x){x$cor_df}) %>%
  dplyr::bind_rows(.id = "region") %>%
  dplyr::filter(feat_A_type != feat_B_type,
                feat_A_type == "RNA",
                feat_B_type == "Metabo")

all_reg_corr_input = lapply(all_reg_corr, function(x){x$cor_input})
names(all_reg_corr_input) = names(all_reg_corr)



# Correlation per cell type region --------------------------------------------------

all_CT_regions = metabo_RNA_seurat@meta.data$CT_Region %>%
  unique() %>% na.omit()
all_CT_regions = all_CT_regions[all_CT_regions != "Uncertain"]

pb = txtProgressBar(min = 0, max = length(all_CT_regions), style = 3)

all_CT_reg_corr = list()
for (i in 1:length(all_CT_regions)){
  all_CT_reg_corr[[all_CT_regions[i]]] = calc_rna_metabo_corr(seurat_obj = metabo_RNA_seurat,
                                                        metadata_filters = list("CT_Region" = all_CT_regions[i]),
                                                        corr_method = "spearman",
                                                        min_nonzero_cells = 30,
                                                        feat_nonzero_prop = 0.5,
                                                        metabo_LOD = 0,
                                                        nonzero_only = F)
  setTxtProgressBar(pb, i)
}

all_CT_reg_corr_df = lapply(all_CT_reg_corr, function(x){x$cor_df}) %>%
  dplyr::bind_rows(.id = "CT_Region") %>%
  dplyr::filter(feat_A_type != feat_B_type,
                feat_A_type == "RNA",
                feat_B_type == "Metabo")

all_CT_reg_corr_input = lapply(all_CT_reg_corr, function(x){x$cor_input})
names(all_CT_reg_corr_input) = names(all_CT_reg_corr)




# Correlation per cell type -----------------------------------------------
all_CT = metabo_RNA_seurat@meta.data$Cell_type_fine %>%
  unique() %>% na.omit() %>% as.character()

pb = txtProgressBar(min = 0, max = length(all_CT), style = 3)

all_CT_corr = list()
for (i in 1:length(all_CT)){
  all_CT_corr[[all_CT[i]]] = calc_rna_metabo_corr(seurat_obj = metabo_RNA_seurat,
                                                        metadata_filters = list("Cell_type_fine" = all_CT[i]),
                                                        corr_method = "spearman",
                                                        min_nonzero_cells = 30,
                                                        feat_nonzero_prop = 0.5,
                                                        metabo_LOD = 0,
                                                        nonzero_only = F)
  setTxtProgressBar(pb, i)
}

all_CT_corr_input = lapply(all_CT_corr, function(x){x$cor_input})
names(all_CT_corr_input) = names(all_CT_corr)

all_CT_corr_df = lapply(all_CT_corr, function(x){x$cor_df}) %>%
  dplyr::bind_rows(.id = "Cell_type") %>%
  dplyr::filter(feat_A_type != feat_B_type,
                feat_A_type == "RNA",
                feat_B_type == "Metabo")

# Correlation per metabotype -----------------------------------------------
all_MT = metabo_RNA_seurat@meta.data$Metabotype %>%
  unique() %>% na.omit() %>% as.character()

pb = txtProgressBar(min = 0, max = length(all_MT), style = 3)

all_MT_corr = list()
for (i in 1:length(all_MT)){
  all_MT_corr[[all_MT[i]]] = calc_rna_metabo_corr(seurat_obj = metabo_RNA_seurat,
                                                  metadata_filters = list("Metabotype" = all_MT[i]),
                                                  corr_method = "spearman",
                                                  min_nonzero_cells = 30,
                                                  feat_nonzero_prop = 0.5,
                                                  metabo_LOD = 0,
                                                  nonzero_only = F)
  setTxtProgressBar(pb, i)
}

all_MT_corr_input = lapply(all_MT_corr, function(x){x$cor_input})
names(all_MT_corr_input) = names(all_MT_corr)

all_MT_corr_df = lapply(all_MT_corr, function(x){x$cor_df}) %>%
  dplyr::bind_rows(.id = "Metabotype") %>%
  dplyr::filter(feat_A_type != feat_B_type,
                feat_A_type == "RNA",
                feat_B_type == "Metabo")


# Save corr results -------------------------------------------------------

saveRDS(list("All_cells" = all_cells_corr_df,
             "Region" = all_reg_corr_df,
             "Celltype" = all_CT_corr_df,
             "Cell_type_and_Region" = all_CT_reg_corr_df,
             "Metabotype" = all_MT_corr_df,
             "Niches" = niches_corr,
             "SG_leiden" = sg_corr_df),
        "../../../Data/Figure 4/M1_S1_S1_Context_specific_corr_results_updated_annots.rds")

saveRDS(list("All_cells" = all_cells_corr_input,
             "Region" = all_reg_corr_input,
             "Celltype" = all_CT_corr_input,
             "Cell_type_and_Region" = all_CT_reg_corr_input,
             "Metabotype" = all_MT_corr_input,
             "Niches" = niches_corr_input,
             "SG_leiden" = sg_corr_input),
        "../../../Data/Figure 4/M1_S1_S1_Context_specific_corr_input_cells_updated_annots.rds")

# Supp Table 17 (complete gene-metabolite correlations) is now generated by
# Analysis_scripts/Downstream/Generate_Supp_Tables.R from the saved Data/ object above.


# Pseudobulk niches correlation -------------------------------------------------------
pseudobulk_seurat = metabo_RNA_seurat
pseudobulk_seurat@graphs = list()
niches_pseudobulk_seurat = AggregateExpression(object = pseudobulk_seurat,
                                   assays = c("Metabo","RNA"), return.seurat = T,
                                   group.by = "Niche")

niches_metabo = as.matrix(GetAssayData(niches_pseudobulk_seurat,
                                       layer = "data",
                                       assay = "Metabo")) %>%
  t() %>%
  as.data.frame()
colnames(niches_metabo) = paste0("metabo_", colnames(niches_metabo))
niches_RNA = as.matrix(GetAssayData(niches_pseudobulk_seurat,
                                    layer = "data",
                                    assay = "RNA")) %>%
  t() %>% as.data.frame()
colnames(niches_RNA) = paste0("rna_", colnames(niches_RNA))

niches_combined_mat = cbind(niches_RNA, niches_metabo) %>% as.matrix()
spearman_cor = cor(niches_combined_mat, method = "spearman") %>% as.data.frame()
pearson_cor = cor(niches_combined_mat, method = "pearson") %>% as.data.frame()
spearman_long = spearman_cor %>%
  rownames_to_column("feat_A") %>%
  pivot_longer(values_to = "coef", names_to = "feat_B", -feat_A) %>%
  dplyr::mutate(feat_A_type = sub("_.*", "", feat_A),
                feat_B_type = sub("_.*", "", feat_B)) %>%
  dplyr::filter(feat_A_type != feat_B_type,
                feat_A_type == "rna",
                feat_B_type == "metabo")

pearson_long = pearson_cor %>%
  rownames_to_column("feat_A") %>%
  pivot_longer(values_to = "coef", names_to = "feat_B", -feat_A) %>%
  dplyr::mutate(feat_A_type = sub("_.*", "", feat_A),
                feat_B_type = sub("_.*", "", feat_B)) %>%
  dplyr::filter(feat_A_type != feat_B_type,
                feat_A_type == "rna",
                feat_B_type == "metabo") %>%
  dplyr::rename("pearson_coef" = "coef")

all_corr = spearman_long %>% dplyr::left_join(pearson_long)

saveRDS(all_corr, "../../../Data/Figure 4/Niches_pseudobulk_corr_results.rds")
saveRDS(niches_combined_mat, "../../../Data/Figure 4/Niches_pseudobulk_corr_input.rds")

corr_wide_mat = all_corr %>%
  dplyr::select(feat_A, feat_B, coef) %>%
  pivot_wider(names_from = feat_A, values_from = coef,values_fill = 0) %>%
  column_to_rownames("feat_B") %>%
  as.matrix()

Heatmap(corr_wide_mat)

