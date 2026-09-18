library(tidyverse)
library(ggpubr)
library(Seurat)
library(grid)
library(ggridges)

source("../Utils/utils.R")
source("../Utils/Plotting.R")
source("../Utils/colors.R")
source("../Utils/Figure_analysis.R")

main_plot_outdir = "../../../../Manuscript/Figures/Main/Figure 4"
supp_plot_outdir = "../../../../Manuscript/Figures/Supp/Figure 4"

# Load data --------------------------------------------------------------------
all_corr_res = readRDS("../Data/Figure 4/M1_S1_S1_Context_specific_corr_results_updated_annots.rds")
metabo_RNA_seurat = readRDS("../Data/Figure 1/updated_Metabo_RNA_combined_seurat.rds")
brain_ontology = read.csv("../Data/Figure 3/Brain_ontology_df.csv")

cor_input = readRDS("../Data/Figure 4/M1_S1_S1_Context_specific_corr_input_cells_updated_annots.rds")

rxn_rxn_dists = readRDS("../Data/Figure 4/gene_metabo_pathway_rxn_dists.rds")

niches_corr = readRDS("../Data/Figure 4/Per_niche_gene_metabo_corr.rds")
niches_pseudobulk_corr = readRDS("../Data/Figure 4/Niches_pseudobulk_corr_results.rds")

# Replace formulas with mol names -----------------------------------------
metabolite_mapping_res = map_metabolite_names(metabo_RNA_seurat)
metabo_RNA_seurat = metabolite_mapping_res$seurat_obj
metabo_molec_mapping = metabolite_mapping_res$mapping

# Plot coef facet region --------------------------------------------------
Thalamas_region_only = brain_ontology$acronym[
  str_detect(brain_ontology$structure_id_path, "549/")
  ]

cor_regions = all_corr_res$Region %>%
  dplyr::mutate(feat_A = sub(".*rna_", "", feat_A),
                feat_B = sub(".*metabo_", "", feat_B)) %>%
  dplyr::filter(abs(coef) > 0.5,
                n_cells_cor >= 30,
                Region %in% Thalamas_region_only) %>%
  dplyr::left_join(metabo_molec_mapping[,c("ion", "mol_name_one_to_one")], by = c("feat_B" = "ion")) %>%
  dplyr::select(-feat_B) %>%
  dplyr::rename("feat_B" = "mol_name_one_to_one") %>%
  dplyr::filter(!is.na(feat_B))

cor_CT_regions = all_corr_res$Cell_type_and_Region %>%
  tidyr::separate(CT_Region, into = c("Cell_type", "Region"), sep = "[.]", remove = F) %>%
  dplyr::mutate(feat_A = sub(".*rna_", "", feat_A),
                feat_B = sub(".*metabo_", "", feat_B)) %>%
  dplyr::filter(abs(coef) > 0.5,
                n_cells_cor >= 30,
                Region %in% Thalamas_region_only) %>%
  dplyr::filter(feat_B %in% metabo_molec_mapping$mol_name_one_to_one)

val_example_coef = cor_CT_regions %>%
  dplyr::filter(Region == "VAL")
val_example_coef$Region = brain_ontology$name[brain_ontology$acronym == "VAL"] %>% unique()
p1 = ggplot(val_example_coef, aes(x = feat_A,
                        y = feat_B,
                        fill = coef)) +
  geom_tile(width = 0.8, height = 0.8) +
  geom_text(aes(label = round(coef, 1)), color = "black") +
  facet_grid(Region ~ Cell_type, scales = "free") +
  scale_fill_gradient2(low = scale_cols$div$low,
                       mid = scale_cols$div$mid,
                       high = scale_cols$div$high, midpoint = 0) +
  theme_pubr() +
  rotate_x_text(angle = 60) +
  theme(strip.text = element_text(size = 12, face = "bold"),
        text = element_text(size = 16)) +
  xlab("") +
  ylab("") +
  labs(fill = "Spearman coefficient")

ggsave(file.path(main_plot_outdir, "VAL_coef_celltypes_example.pdf"),
       p1, width = 8, height = 7)


# Plot coef distr comparison ----------------------------------------------

corr_res_subset = all_corr_res[which(!names(all_corr_res) %in% c("All_cells", "Cell_type_and_Region"))]
same_size_coef = lapply(corr_res_subset, function(x){
  x %>%
    dplyr::filter(feat_B %in% paste0("metabo_",metabo_molec_mapping$ion)) %>%
    dplyr::group_by(feat_A, feat_B) %>%
    dplyr::slice_max(abs(coef), n = 1, with_ties = F) %>%
    dplyr::select(coef)
}) %>%
  dplyr::bind_rows(.id = "Category")

same_size_coef_global = lapply(list("All Cells" = all_corr_res$All_cells,
                                    "All Niches" = niches_pseudobulk_corr), function(x){
  x %>%
    dplyr::filter(feat_B %in% paste0("metabo_",metabo_molec_mapping$ion)) %>%
    dplyr::group_by(feat_A, feat_B) %>%
    dplyr::slice_max(abs(coef), n = 1, with_ties = F) %>%
    dplyr::select(coef)
}) %>%
  dplyr::bind_rows(.id = "Category")

# Plot comparative density plot
p1 <- ggplot(same_size_coef, aes(x = coef,
                                 y = Category,
                                 fill = Category,
                                 color = Category)) +
  geom_density_ridges(alpha = 0.2, size = 1, scale = 1) +
  theme_minimal() +
  labs(x = "Spearman Coefficient",
       y = "",
       fill = "Category",
       color = "Category") +
  theme_pubr()
p1 = p1 +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  scale_color_manual(values = assign_colors(unique(same_size_coef$Category),
                                            cols = my_cols$corr_labels_histo)) +
  scale_fill_manual(values = assign_colors(unique(same_size_coef$Category),
                                           cols = my_cols$corr_labels_histo)) +
  theme(text = element_text(size = 16))

ggsave(file.path(main_plot_outdir, "Corr_distribution_per_context.pdf"),
       p1, width = 7, height = 4)

p2 <- ggplot(same_size_coef_global, aes(x = coef,
                                        y = Category,
                                        fill = Category,
                                        color = Category)) +
  geom_density_ridges(alpha = 0.2, size = 1, scale = 1) +
  theme_minimal() +
  labs(x = "Spearman Coefficient",
       y = "",
       fill = "Category",
       color = "Category") +
  theme_pubr()
p2 = p2 %>% clean_labels() +
  scale_color_manual(values = assign_colors(unique(same_size_coef_global$Category),
                                            cols = my_cols$highlight_colors)) +
  scale_fill_manual(values = assign_colors(unique(same_size_coef_global$Category),
                                           cols = my_cols$highlight_colors)) +
  theme(text = element_text(size = 16))

ggsave(file.path(main_plot_outdir, "Corr_distribution_global.pdf"),
       p2, width = 6, height = 4)


# Plot spatial panels with correlation ------------------------------------

p = plot_spatial_cor_res(seurat_obj = metabo_RNA_seurat,
                         Region_name = "VAL",
                         cor_input_data = cor_input,
                         rna_feat = "rna_Hk1",
                         metabo_feat = "metabo_D-Glucose*+Cl",
                         cell_type = "Oligodendrocytes",
                         coef = -0.581709518,
                         convert_sf_to_name = F)

p2 = ggarrange(plotlist = p[-4],nrow = 2, ncol = 2)
export_plot(plot = p2 %>% clean_labels(),
            path = file.path(supp_plot_outdir, "VAL_corr_spatial_panels_example.pdf"),
            sideways = T)


# Overlap with reaction-reaction network dists ----------------------------
sf_name_unique = read.csv("../Data/General/updated_annots_w_1_1_mapping.csv") %>%
  dplyr::select(ion, mol_name_one_to_one) %>%
  dplyr::mutate(ion = sub("[+-].*", "", ion)) %>%
  dplyr::distinct()


tmp = tempfile()
utils::download.file(
  "https://rest.kegg.jp/list/pathway/mmu",
  destfile = tmp)
kegg_mouse_pathways = read.delim(tmp, header = F)
colnames(kegg_mouse_pathways) = c("kegg_pathway_id", "kegg_pathway_name")

all_dists = rxn_rxn_dists %>% dplyr::filter(!is.na(dist))

pathway_mats = build_gene_metabo_pathway_mats(niches_corr, sf_name_unique, all_dists, kegg_mouse_pathways)
cor_info_mat = pathway_mats$cor_info_mat
pathway_info_mat = pathway_mats$pathway_info_mat
pathway_info_label_mat = pathway_mats$pathway_info_label_mat


cor_heatmap = ComplexHeatmap::Heatmap(cor_info_mat,cluster_rows = F,
                            cluster_columns = F,
                            column_title_gp = gpar(fontsize = 15),
                            col = circlize::colorRamp2(c(-1, 0, 1), c(scale_cols$div$low,
                                                            scale_cols$div$mid,
                                                            scale_cols$div$high)),
                            column_names_side = "top",
                            column_names_rot = 45,
                            cell_fun = function(j, i, x, y, width, height, fill) {
                              if(abs(cor_info_mat[i, j]) > 0){
                                grid.text(sprintf("%.1f", cor_info_mat[i, j]), x, y, gp = gpar(fontsize = 12))
                              }
                              else{
                                grid.rect(x,y,width,height,gp = gpar(fill = "white"))
                              }
                            },name = "Niche correlation",
                            heatmap_legend_param = list(title_gp = gpar(fontsize = 15)))
pathway_heatmap = ComplexHeatmap::Heatmap(pathway_info_mat,cluster_rows = F,
                            cluster_columns = F,
                            col = c("white", "black"),
                            column_names_rot = 45,
                            row_title_gp = gpar(fontsize = 15),
                            column_title_gp = gpar(fontsize = 15),
                            cell_fun = function(j, i, x, y, width, height, fill) {
                              if(pathway_info_mat[i, j] > 0){
                                grid.text(sprintf(pathway_info_label_mat[i, j]), x, y, gp = gpar(fontsize = 10,
                                                                                                 col = "white"))
                              }
                              else{
                                grid.rect(x,y,width,height,gp = gpar(fill = "white"))
                              }
                            },
                            name = "Pathway information",
                            heatmap_legend_param = list(title_gp = gpar(fontsize = 15)))

final_plot = cor_heatmap + pathway_heatmap

# ggsave(file.path(main_plot_outdir, "Metabo_gene_corr_pathway_KEGG_heatmaps.pdf"),
#        final_plot, width = 7, height = 4)
#This plot was saved manually as 7 by 13 landscape because it was big
# export_plot(plot = final_plot,
#             path = file.path(main_plot_outdir, "Metabo_gene_corr_pathway_KEGG_heatmaps.pdf"),
#             sideways = T)


niche_plot_seurat = metabo_RNA_seurat
Idents(niche_plot_seurat) = "Niche"

sel_niches = pathway_mats$niches_of_interest
sel_cells = vector("list", length(sel_niches)) %>% setNames(sel_niches)
for (i in sel_niches){
  sel_cells[[i]] = WhichCells(niche_plot_seurat, idents = i)
}

sel_spatial_niches_plot = SpatialDimPlot(object = niche_plot_seurat,
               cells.highlight = sel_cells,
               cols.highlight = c(c(rep(my_cols$highlight_colors[2],7),
                                    my_cols$highlight_colors[1],
                                    rep(my_cols$highlight_colors[2],5)),"lightgrey"),
               stroke = NA,
               label = T,
               pt.size.factor = 1,
               label.size = 4,
               label.box = F,
               label.color = "black") +
  theme(legend.position = "none") +
  coord_fixed()

ggsave(file.path(main_plot_outdir, "Selected_Niches_spatial_plot_for_corr_pathway_heatmap.pdf"),
       sel_spatial_niches_plot, width = 6, height = 6)

# Plot pred accuracy RMSE -------------------------------------------------
# Matches Run_MISTY_gene_metabo.R's output (re-run that script if these are missing)
with_intra_dir = "../Data/Figure 4/updated_with_intra"
without_intra_dir = "../Data/Figure 4/updated_without_intra"

RMSE_final = misty_rmse_comparison_plots(with_intra_dir, without_intra_dir)
export_plot(plot = RMSE_final,
            path = file.path(supp_plot_outdir, "RMSE_misty_plots_and_contribution.pdf"),
            sideways = T)

# Manuscript-text numbers (figure_4_numbers.csv, figure_4_VAL_example_n_cells_per_celltype.csv)
# are generated by Fig_scripts/Generate_Figure_Numbers.R.
# Note: the Aldh2-specific example correlations (alde_metabo1_rho / alde_metabo2_rho) referenced
# in the manuscript's Aldh2 pathway illustration (Fig 4F) are not computed by this script -- that
# panel is not generated here, flag for manual verification.
