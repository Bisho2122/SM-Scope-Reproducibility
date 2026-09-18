library(S2IsoMEr)
library(tidyverse)

source("../../Utils/utils.R")

# Load Data ---------------------------------------------------------------
DE_res = readRDS("../../Data/Figure 2/updated_sc_DE_markers.rds")
all_misty_res = readRDS("../../Data/Figure 3/updated_MISTY_metabo_region_result.rds")
metabo_RNA_seurat = readRDS("../../Data/Figure 1/updated_Metabo_RNA_combined_seurat.rds")

metabo_subclass_bg = Load_background(mol_type = "Metabo",
                                     bg_type = "sub_class",
                                     feature_type = "sf")

metabo_pathway_bg = Load_background(mol_type = "Metabo",
                                    bg_type = "pathway",
                                    feature_type = "sf")

# Replace formulas with mol names -----------------------------------------
metabo_molec_mapping = read.csv("../../Data/General/updated_annots_w_1_1_mapping.csv")
metabo_molec_mapping = metabo_molec_mapping[metabo_molec_mapping$ion %in% rownames(metabo_RNA_seurat@assays$Metabo),]

metabo_molec_mapping$adduct = str_extract(metabo_molec_mapping$ion,
                                          "(\\+|\\-).*")
metabo_molec_mapping$mol_name_one_to_one = paste0(metabo_molec_mapping$mol_name_one_to_one, metabo_molec_mapping$adduct)

DE_res$Metabo[c("Celltypes", "Metabotype","Region", "CT_Region")] <- lapply(
  DE_res$Metabo[c("Celltypes", "Metabotype","Region", "CT_Region")],
  map_sf_to_names
)


# Get markers per metabotype ----------------------------------------------
to_remove_metabo = c("prontosil-H")
DE_metabotypes = DE_res$Metabo$Metabotype %>%
  dplyr::filter(avg_log2FC >= 1,
                p_val_adj < 0.05,
                pct.1 >= 0.5,
                ! gene %in% to_remove_metabo)

all_metabotypes = unique(DE_metabotypes$cluster) %>% as.character()

metabotypes_markers = list()
for (r in all_metabotypes){
  de_mols = DE_metabotypes$gene[DE_metabotypes$cluster == r]
  de_sf = metabo_molec_mapping$ion[metabo_molec_mapping$mol_name_one_to_one %in% de_mols]
  de_sf = sub("[+-].*", "", de_sf)
  metabotypes_markers[[r]] = de_sf %>% unique()
}

# Get markers per region --------------------------------------------------
to_remove_metabo = c("prontosil-H")
DE_regs = DE_res$Metabo$Region %>%
  dplyr::filter(avg_log2FC >= 1,
                p_val_adj < 0.05,
                pct.1 >= 0.5,
                cluster != "Uncertain",
                ! gene %in% to_remove_metabo)

all_regs = unique(DE_regs$cluster) %>% as.character()

region_markers = list()
for (r in all_regs){
  de_mols = DE_regs$gene[DE_regs$cluster == r]
  de_sf = metabo_molec_mapping$ion[metabo_molec_mapping$mol_name_one_to_one %in% de_mols]
  de_sf = sub("[+-].*", "", de_sf)
  region_markers[[r]] = de_sf %>% unique()
}

# Prepare MISTY region MSEA input -----------------------------------------
misty_regions_scores = all_misty_res$scenario_3$Regions$importances.aggregated %>%
  dplyr::filter(view == "Metabo") %>%
  dplyr::mutate(pred_sf = sub("_.*", "", Predictor)) %>%
  dplyr::group_by(Target, pred_sf) %>%
  dplyr::slice_max(Importance, n = 1, with_ties = F)

all_regs = unique(misty_regions_scores$Target)
misty_msea_input = list()
for (r in all_regs){
  filt = misty_regions_scores[misty_regions_scores$Target == r,]
  vals = filt$Importance
  names(vals) = filt$pred_sf
  vals = vals[order(vals)]
  misty_msea_input[[r]] = vals
}
# Run and plot ORA Region -----------------------------------------------------------------
# univ = unique(unlist(region_markers))
univ = sub("[+-].*", "", metabo_molec_mapping$ion) %>% unique()

subclass_ORA_res = Run_simple_ORA(marker_list = region_markers,
                                  background = metabo_subclass_bg,
                                  custom_universe = univ,
                                  alpha_cutoff = 1,
                                  min_intersection = 2)

# pathway_ORA_res = Run_simple_ORA(marker_list = region_markers,
#                                   background = metabo_pathway_bg,
#                                   custom_universe = univ,
#                                   alpha_cutoff = 1,
#                                   min_intersection = 2)


subclass_ORA_res$q.value = subclass_ORA_res$p_value
# pathway_ORA_res$q.value = pathway_ORA_res$p_value

subclass_ORA_res = subclass_ORA_res %>%
  dplyr::filter(q.value <= 0.05)

# disease_terms = c()
# for (i in c("syndrome", "deficiency", "disease",
#             "Deficiency", "Diseases", "Disease",
#             "Syndrome", "disorder", "disorders",
#             "Action Pathway")){
#   disease_terms = c(disease_terms,
#                     str_subset(unique(pathway_ORA_res$Term),
#                                i))
# }

# pathway_ORA_res = pathway_ORA_res %>%
#   dplyr::filter(q.value <= 0.05,
#                 !Term %in% disease_terms)

ORA_region = subclass_ORA_res


p = dotplot_ORA(ORA_res = subclass_ORA_res) +
  ggpubr::theme_pubr(legend = "right") +
  ggpubr::rotate_x_text(45)

# Run and plot ORA metabotypes -----------------------------------------------------------------
# univ = unique(unlist(region_markers))
univ = sub("[+-].*", "", metabo_molec_mapping$ion) %>% unique()

subclass_ORA_res = Run_simple_ORA(marker_list = metabotypes_markers,
                                  background = metabo_subclass_bg,
                                  custom_universe = univ,
                                  alpha_cutoff = 1,
                                  min_intersection = 2)

subclass_ORA_res$q.value = subclass_ORA_res$p_value

subclass_ORA_res = subclass_ORA_res %>%
  dplyr::filter(q.value <= 0.05)

ORA_metabotypes = subclass_ORA_res


p1 = dotplot_ORA(ORA_res = subclass_ORA_res) +
  ggpubr::theme_pubr(legend = "right") +
  ggpubr::rotate_x_text(45)


# Save ORA results --------------------------------------------------------
saveRDS(ORA_metabotypes, "../../Data/Figure 2/updated_ORA_metabotypes.rds")
saveRDS(ORA_region, "../../Data/Figure 3/updated_ORA_regions.rds")

# Supp Tables 11 and 14 (metabotype/region ORA) are now generated by
# Analysis_scripts/Downstream/Generate_Supp_Tables.R from these saved Data/ objects.


# Run MSEA MISTY Regions --------------------------------------------------
# NOTE: this analysis's results are not included in the final manuscript. Kept in place
# (not archived) per the author's choice, but flagged here since nothing else in the repo
# consumes MSEA_misty_regions_results.rds.
library(fgsea)
metabo_subclass_bg = Load_background(mol_type = "Metabo",
                                     bg_type = "sub_class",
                                     feature_type = "sf")
# msea_subclass_bg = lapply(metabo_subclass_bg, function(x){
#   x[x %in% misty_regions_scores$pred_sf]
# })
# msea_subclass_bg = msea_subclass_bg[lengths(msea_subclass_bg) > 0]

msea_misty_res = list()
for (i in 1:length(misty_msea_input)){
  reg = names(misty_msea_input)[i]
  msea_res = fgsea::fgseaMultilevel(pathways = metabo_subclass_bg,
                                    stats = misty_msea_input[[i]],
                                    minSize = 3)
  msea_misty_res[[reg]] = msea_res
}
saveRDS(msea_misty_res, "../../Data/Figure 3/MSEA_misty_regions_results.rds")

msea_misty_res = msea_misty_res %>% dplyr::bind_rows(.id = "region")
msea_misty_res = msea_misty_res %>%
  dplyr::filter(pval < 0.05) %>%
  dplyr::arrange(NES)

mat <- msea_misty_res %>%
  select(pathway, region, NES) %>%
  pivot_wider(names_from = region, values_from = NES, values_fill = 0) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

# Step 2: Hierarchical clustering
row_clust <- hclust(dist(mat))       # for pathways
col_clust <- hclust(dist(t(mat)))    # for groups

# Step 3: Use factor levels to preserve clustering order
pathway_levels <- rownames(mat)[row_clust$order]
group_levels <- colnames(mat)[col_clust$order]

plot_data <- msea_misty_res %>%
  filter(pathway %in% pathway_levels, region %in% group_levels) %>%
  mutate(
    pathway = factor(pathway, levels = pathway_levels),
    region = factor(region, levels = group_levels)
  )

ggplot(plot_data, aes(x = region, y = pathway)) +
  geom_point(aes(color = NES, size = -log10(pval))) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "NES") +
  scale_size(range = c(2, 6), name = "-log10(pvalue)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) +
  ylab("Pathway") + xlab("Region") +
  ggtitle("FGSEA Dot Plot: NES by Region and Pathway")
