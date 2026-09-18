# Generates every manuscript-text number CSV in Output/Numbers/ from a single place, instead of
# each Figure_N.R computing and writing its own numbers inline. Mirrors
# Analysis_scripts/Downstream/Generate_Supp_Tables.R's structure (one section per figure) and
# reuses the same analysis functions each Figure_N.R calls for its plots (Utils/Figure_analysis.R),
# so the values here can't drift from what the figures themselves show.

library(tidyverse)
library(Seurat)
require(DIALOGUE)

source("../Utils/utils.R")
source("../Utils/Plotting.R")
source("../Utils/colors.R")
source("../Utils/Figure_analysis.R")

out_dir = "../Output/Numbers"

# Figure 2 -----------------------------------------------------------------
metabo_RNA_seurat = readRDS("../Data/Figure 1/updated_Metabo_RNA_combined_seurat.rds")
DE_res = readRDS("../Data/Figure 2/updated_sc_DE_markers.rds")

figure_2_numbers = compute_figure2_numbers(metabo_RNA_seurat, DE_res)
write.csv(figure_2_numbers, file.path(out_dir, "figure_2_numbers.csv"), row.names = FALSE)

# Figure 3 -----------------------------------------------------------------
SNR_fisher_res = readRDS("../Data/Figure 3/SNR_region_fisher_res.rds")

figure_3_numbers = compute_figure3_numbers(metabo_RNA_seurat, SNR_fisher_res)
write.csv(figure_3_numbers, file.path(out_dir, "figure_3_numbers.csv"), row.names = FALSE)

# Figure 4 -----------------------------------------------------------------
all_corr_res = readRDS("../Data/Figure 4/M1_S1_S1_Context_specific_corr_results_updated_annots.rds")
rxn_rxn_dists = readRDS("../Data/Figure 4/gene_metabo_pathway_rxn_dists.rds")
brain_ontology = read.csv("../Data/Figure 3/Brain_ontology_df.csv")

metabolite_mapping_res = map_metabolite_names(metabo_RNA_seurat, rename_assay_rownames = FALSE)
metabo_molec_mapping = metabolite_mapping_res$mapping

all_dists = rxn_rxn_dists %>% dplyr::filter(!is.na(dist))

figure_4_res = compute_figure4_numbers(all_corr_res, all_dists, brain_ontology, metabo_molec_mapping)
write.csv(figure_4_res$numbers, file.path(out_dir, "figure_4_numbers.csv"), row.names = FALSE)
write.csv(figure_4_res$n_cells_in_VAL_example,
          file.path(out_dir, "figure_4_VAL_example_n_cells_per_celltype.csv"), row.names = FALSE)

# Figure 5 -----------------------------------------------------------------
dialogue_res = readRDS("../Data/Figure 5/updated_DIALOGUE_res/updated_CT_5_main_out.rds")
all_DE_res = readRDS("../Data/Figure 2/updated_sc_DE_markers.rds")
spatial_niches = readRDS("../Data/Figure 1/spatial_niches_regions.rds")

figure_5_numbers = compute_figure5_numbers(metabo_RNA_seurat, dialogue_res, all_DE_res, spatial_niches)
write.csv(figure_5_numbers, file.path(out_dir, "figure_5_numbers.csv"), row.names = FALSE)

# Figure 6 -----------------------------------------------------------------
all_sg_res = readRDS("../Data/Figure 6/updated_SG_eval_optim_res.rds")

sel_sg_res = all_sg_res
sel_sg_res$SG_clusts_metadata$Cell_type_fine = str_replace_all(sel_sg_res$SG_clusts_metadata$Cell_type_fine,
                                                               c("Excit_CTX" = "Excitatory_VGLUT1",
                                                                 "Inhib_CTX" = "Inhibitory_cortex",
                                                                 "Excit_VGLUT2" = "Excitatory_VGLUT2",
                                                                 "Excit_HPF" = "Excitatory_HPF",
                                                                 "Inhib_Interneurons" = "Inhibitory_interneurons",
                                                                 "Ventricular systems Glial cells" = "Ependymal cells"))
sel_sg_res$SG_clusts_metadata$Cell_type_coarse = str_replace_all(sel_sg_res$SG_clusts_metadata$Cell_type_coarse,
                                                               c("Excit" = "Excitatory",
                                                                 "Inhib" = "Inhibitory",
                                                                 "Ventricular systems Glial cells" = "Ependymal cells"))

figure_6_res = compute_figure6_numbers(metabo_RNA_seurat, sel_sg_res)
write.csv(figure_6_res$numbers, file.path(out_dir, "figure_6_numbers.csv"), row.names = FALSE)
write.csv(figure_6_res$dominant_SG_per_region,
          file.path(out_dir, "figure_6_dominant_SG_per_region.csv"), row.names = FALSE)
