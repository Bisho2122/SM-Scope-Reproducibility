# Generates every supplementary table CSV in Output/Supp_Tables/ from a single place. Reads
# existing Data/ objects and reformats them -- doesn't recompute analyses, except the
# region-parent metabolite markers table (needs a FindAllMarkers() re-run, see the "not
# referenced in the manuscript" section at the end). Table numbers match the manuscript's
# supplementary table order. Tables 16-18 are sourced from
# Analysis_scripts/Upstream/Resolve_Probe_Panel_Design/Design_Resolve_probe_panel.R.

library(tidyverse)
library(Seurat)

source("../../Utils/utils.R")

# Set to your own cluster storage path -- needed for Table 1, Tables 16-17, and the Baysor
# cell-stats table only, everything else reads from Data/.
cluster_base_dir = Sys.getenv("SPATIAL_OMICS_CLUSTER_DIR", "/path/to/your/cluster/storage/Spatial_omics")
out_dir = "../../Output/Supp_Tables"

# Shared objects, loaded once --------------------------------------------------------------
metabo_RNA_seurat = readRDS("../../Data/Figure 1/updated_Metabo_RNA_combined_seurat.rds")

metabo_molec_mapping = read.csv("../../Data/General/updated_annots_w_1_1_mapping.csv")
metabo_molec_mapping$adduct = str_extract(metabo_molec_mapping$ion, "(\\+|\\-).*")
metabo_molec_mapping$mol_name_one_to_one = paste0(metabo_molec_mapping$mol_name_one_to_one, metabo_molec_mapping$adduct)

DE_res = readRDS("../../Data/Figure 2/updated_sc_DE_markers.rds")
brain_ontology = read.csv("../../Data/Figure 3/Brain_ontology_df.csv")

# Loaded here (rather than in Table 6's/13's own section) because each is also needed by a
# table that now sits earlier in the manuscript-ordered script (Table 6 needs it before
# Table 19 generates it; Table 13 needs it before Table 20 generates it).
final_regions = readRDS("../../Data/Figure 1/Final_regions_atlas.rds")
misty_niches_res = readRDS("../../Data/Figure 4/updated_MISTY_Niches_pseudobulk_gene_metabo.rds")

sig_filter = function(df) {
  df %>% dplyr::filter(abs(avg_log2FC) > 1, p_val_adj < 0.05)
}


# Table 1: METASPACE annotation results per database, with off-sample QC ---------------------
raw_annots = read.csv(file.path(cluster_base_dir, "Results/Upstream/METASPACE_annotations/2022-01-11_23h33m36s.csv"))

qc_dir = file.path(cluster_base_dir, "Results/Upstream/METASPACE_annotations/Annotation_QC")
qc_files = list.files(qc_dir, pattern = "^Annot_QC_2022.*\\.csv$")
db_name_fix = c("RAMP_DB" = "RaMP_db_brain")

off_sample_qc = qc_files %>%
  purrr::map_dfr(~ read.csv(file.path(qc_dir, .x))) %>%
  dplyr::mutate(DB_name = dplyr::recode(DB_name, !!! db_name_fix))

supp_table_1 = raw_annots %>%
  dplyr::filter(database != "HMDB-endogenous") %>%
  dplyr::left_join(off_sample_qc, by = c("sumFormula" = "formula", "adduct" = "adduct", "database" = "DB_name")) %>%
  dplyr::select(database, sumFormula, adduct, ion, mz, msmScore, rhoSpatial, rhoSpectral, rhoChaos,
                fdrLevel, off_sample_label, off_sample_prob)

write.csv(supp_table_1, file.path(out_dir, "Supp_Table_1_METASPACE_annotation_results.csv"), row.names = FALSE)


# Table 2: Cell-type annotation details (Leiden -> fine type, markers) ----------------------
seurat_res = readRDS("../../Data/Figure 1/Final_seurat_transc_updated.rds")$seurat_res

# Same cell-type renaming as Create_multi_assay_seurat.R, applied here directly (rather than
# joining in metabo_RNA_seurat's already-renamed Cell_type_fine) so every Leiden cluster in
# seurat_res$analysis_markers gets a label, including any cluster whose cells didn't all
# survive that script's RNA/Metabo common-cell subsetting.
leiden_to_fine = data.frame(
  seurat_clusters = as.character(seurat_res$Seurat_obj$seurat_clusters),
  Cell_type_fine = as.character(Idents(seurat_res$Seurat_obj))
) %>%
  dplyr::mutate(Cell_type_fine = str_replace_all(Cell_type_fine,
                                                  c("Excit_CTX" = "Excitatory_VGLUT1",
                                                    "Inhib_CTX" = "Inhibitory_cortex",
                                                    "Excit_VGLUT2" = "Excitatory_VGLUT2",
                                                    "Excit_HPF" = "Excitatory_HPF",
                                                    "Inhib_Interneurons" = "Inhibitory_interneurons",
                                                    "Ventricular systems Glial cells" = "Ependymal cells"))) %>%
  dplyr::mutate(Cell_type_coarse = sub("_.*", "", Cell_type_fine)) %>%
  dplyr::distinct()

supp_table_2 = seurat_res$analysis_markers %>%
  dplyr::left_join(leiden_to_fine, by = c("cluster" = "seurat_clusters"))

write.csv(supp_table_2, file.path(out_dir, "Supp_Table_2_celltype_annotation_details.csv"), row.names = FALSE)


# Table 3: Final high-confidence metabolite panel with 1:1 mapping -------------------------
final_ions = rownames(metabo_RNA_seurat@assays$Metabo)
supp_table_3 = metabo_molec_mapping %>%
  dplyr::filter(ion %in% final_ions)

write.csv(supp_table_3, file.path(out_dir, "Supp_Table_3_final_metabolite_panel.csv"), row.names = FALSE)


# Table 4: Metabotype-specific metabolite DE markers ----------------------------------------
supp_table_4 = sig_filter(DE_res$Metabo$Metabotype)

write.csv(supp_table_4, file.path(out_dir, "Supp_Table_4_metabotype_metabolite_markers.csv"), row.names = FALSE)


# Table 5: Metabolic class enrichment (ORA) per metabotype ----------------------------------
ORA_metabotypes = readRDS("../../Data/Figure 2/updated_ORA_metabotypes.rds")

write.csv(ORA_metabotypes, file.path(out_dir, "Supp_Table_5_metabotype_ORA.csv"), row.names = FALSE)


# Table 6: Allen Brain Atlas structures used in this study ----------------------------------
regions_used = unique(final_regions$final_gpd$region)

supp_table_6 = brain_ontology %>%
  dplyr::filter(acronym %in% regions_used) %>%
  dplyr::select(acronym, id, name, parent_structure_id, structure_id_path, depth)

write.csv(supp_table_6, file.path(out_dir, "Supp_Table_6_allen_atlas_structures_used.csv"), row.names = FALSE)


# Table 7: Region-specific metabolite markers (fine-grained) --------------------------------
supp_table_7 = dplyr::bind_rows(list("Region" = sig_filter(DE_res$Metabo$Region),
                                      "CT_Region" = sig_filter(DE_res$Metabo$CT_Region)),
                                 .id = "Grouping") %>%
  dplyr::filter(Grouping == "Region") %>%
  dplyr::select(-Grouping)

write.csv(supp_table_7, file.path(out_dir, "Supp_Table_7_region_fine_metabolite_markers.csv"), row.names = FALSE)


# Table 8: Metabolic class enrichment (ORA) per region ---------------------------------------
ORA_region = readRDS("../../Data/Figure 3/updated_ORA_regions.rds")

write.csv(ORA_region, file.path(out_dir, "Supp_Table_8_region_ORA.csv"), row.names = FALSE)


# Table 9: MISTY region-predictive metabolites ------------------------------------------------
all_misty_res = readRDS("../../Data/Figure 3/updated_MISTY_metabo_region_result.rds")

supp_table_9 = Get_interaction_importance(misty.results = all_misty_res$scenario_3$Regions,
                                           view = "Metabo", cutoff = 2,
                                           trim.measure = "gain.RMSE", trim = 10, clean = T) %>%
  dplyr::select(-nsamples)

write.csv(supp_table_9, file.path(out_dir, "Supp_Table_9_MISTY_region_predictive_metabolites.csv"), row.names = FALSE)


# Table 10: MISTY SNR enrichment ---------------------------------------------------------------
fisher_res_all = readRDS("../../Data/Figure 3/SNR_region_fisher_res.rds")

write.csv(fisher_res_all, file.path(out_dir, "Supp_Table_10_MISTY_SNR_enrichment.csv"), row.names = FALSE)


# Table 11: Gene-metabolite correlations across all contexts ---------------------------------------
all_corr_res = readRDS("../../Data/Figure 4/M1_S1_S1_Context_specific_corr_results_updated_annots.rds")

supp_table_11 = dplyr::bind_rows(all_corr_res, .id = "Context") %>%
  dplyr::filter(n_cells_cor >= 30, abs(coef) >= 0.5)

write.csv(supp_table_11, file.path(out_dir, "Supp_Table_11_gene_metabo_correlations.csv"), row.names = FALSE)


# Table 12: Pathway-contextualized gene-metabolite pairs ---------------------------------------
all_dists = readRDS("../../Data/Figure 4/gene_metabo_pathway_rxn_dists.rds") %>%
  dplyr::filter(!is.na(dist))

write.csv(all_dists, file.path(out_dir, "Supp_Table_12_pathway_gene_metabo_distances.csv"), row.names = FALSE)


# Table 13: MISTY vs correlation comparison -----------------------------------------------------
niches_pb_corr = readRDS("../../Data/Figure 4/Niches_pseudobulk_corr_results.rds")

niches_misty_pb_clean = Get_interaction_importance(misty.results = misty_niches_res,
                                                    view = "gene", cutoff = 1,
                                                    trim.measure = "gain.R2", trim = 1, clean = T)
ion_adduct_mapping = data.frame(ion = rownames(metabo_RNA_seurat@assays$Metabo),
                                ion_no_adduct = sub("[+-]", "_", rownames(metabo_RNA_seurat@assays$Metabo)))

niches_corr_p13 = niches_pb_corr %>%
  dplyr::mutate(feat_A = sub(".*rna_", "", feat_A),
                feat_B = sub(".*metabo_", "", feat_B)) %>%
  dplyr::left_join(ion_adduct_mapping, by = c("feat_B" = "ion")) %>%
  dplyr::select(-feat_B) %>%
  dplyr::rename("feat_B" = "ion_no_adduct")

niches_misty_corr = niches_misty_pb_clean %>%
  dplyr::left_join(niches_corr_p13, by = c("Target" = "feat_B", "Predictor" = "feat_A")) %>%
  dplyr::filter(!is.na(coef))

supp_table_13 = niches_misty_corr %>%
  dplyr::group_by(Target) %>%
  dplyr::slice_max(Importance, n = 3, with_ties = T) %>%
  dplyr::filter(between(coef, 0, 0.299), between(pearson_coef, 0, 0.299), Importance >= 2) %>%
  dplyr::ungroup()

write.csv(supp_table_13, file.path(out_dir, "Supp_Table_13_MISTY_vs_correlation.csv"), row.names = FALSE)


# Table 14: Spatial niches per region -----------------------------------------------------------
spatial_niches = readRDS("../../Data/Figure 1/spatial_niches_regions.rds")
niche_size = spatial_niches %>%
  dplyr::group_by(Niche) %>%
  dplyr::summarise(n = n())
cell_coords = metabo_RNA_seurat@images$image@coordinates %>%
  dplyr::select(x, y) %>%
  rownames_to_column("Cell")

supp_table_14 = spatial_niches %>%
  dplyr::left_join(niche_size, by = "Niche") %>%
  dplyr::left_join(cell_coords, by = "Cell") %>%
  dplyr::rename(n_cells_in_niche = n) %>%
  dplyr::select(-x, -y) %>%
  dplyr::distinct()

write.csv(supp_table_14, file.path(out_dir, "Supp_Table_14_spatial_niches.csv"), row.names = FALSE)


# Table 15: DIALOGUE MCP results + metabolic signatures ------------------------------------------
Res = readRDS("../../Data/Figure 5/CT_5_run_main_out.rds")

supp_table_15 = MCP_signature_to_df(sig_list = Res$R$sig2, min_celltype = 5) %>%
  dplyr::distinct()

write.csv(supp_table_15, file.path(out_dir, "Supp_Table_15_MCP_metabolic_signatures.csv"), row.names = FALSE)

# Table 16: KEGG metabolic pathway list used for probe-panel gene selection --------------------
KEGG_metabolic_pathways = readRDS(file.path(cluster_base_dir, "Data/Resolve_Probe_Panel_Design/KEGG_metabolic_pathways_clean.rds"))

supp_table_16 = KEGG_metabolic_pathways

write.csv(supp_table_16, file.path(out_dir, "Supp_Table_16_KEGG_metabolic_pathways.csv"), row.names = FALSE)


# Table 17: Resolve mouse brain probe panel (reference gene list) ------------------------------
Resolve_probes = readxl::read_xlsx(file.path(cluster_base_dir, "Data/Resolve_Probe_Panel_Design/Resolve_mousebrain_probe_genes.xlsx"))

supp_table_17 = Resolve_probes

write.csv(supp_table_17, file.path(out_dir, "Supp_Table_17_Resolve_probe_panel.csv"), row.names = FALSE)


# Table 18: Final 100-gene Resolve probe panel selection ----------------------------------------
resolve_final_100 = readxl::read_xlsx(file.path(cluster_base_dir, "Data/Resolve_Probe_Panel_Design/final_100.xlsx"))

supp_table_18 = resolve_final_100

write.csv(supp_table_18, file.path(out_dir, "Supp_Table_18_final_100_gene_panel.csv"), row.names = FALSE)


# Table 19: Cell-to-region assignment matrix -------------------------------------------------
supp_table_19 = final_regions$final_gpd %>%
  dplyr::select(Cell = X__index, region, cell_fraction)

write.csv(supp_table_19, file.path(out_dir, "Supp_Table_19_cell_to_region_assignment.csv"), row.names = FALSE)


# Table 20: MISTY gene-to-metabolite prediction ------------------------------------------------
supp_table_20 = Get_interaction_importance(misty.results = misty_niches_res,
                                           view = "gene", cutoff = 1.5,
                                           trim.measure = "gain.R2", trim = 1, clean = T) %>%
  dplyr::select(-nsamples)

write.csv(supp_table_20, file.path(out_dir, "Supp_Table_20_MISTY_gene_metabolite_predictions.csv"), row.names = FALSE)


# Table 21: MCP-metabotype overlap ------------------------------------------------------------
dialogue_res_overlap = readRDS("../../Data/Figure 5/updated_DIALOGUE_res/updated_CT_5_main_out.rds")

MCP_cell_embed = dialogue_res_overlap$full_RA_res$scores %>%
  dplyr::bind_rows() %>%
  rownames_to_column("Cell")
metabotype_meta = metabo_RNA_seurat@meta.data %>%
  rownames_to_column("Cell") %>%
  dplyr::select(Cell, Metabotype)
conf_mat_MCP_metabotype = MCP_cell_embed %>%
  dplyr::left_join(metabotype_meta)
conf_mat_MCP_metabotype = conf_mat_MCP_metabotype[, c("Cell", paste0("MCP", 1:10), "Metabotype")]
conf_mat_MCP_metabotype = conf_mat_MCP_metabotype %>%
  gather(key = "MCP", value = "MCP_score", -Cell, -Metabotype)

MCP_cells_assign = conf_mat_MCP_metabotype %>%
  dplyr::select(Cell, MCP, MCP_score) %>%
  dplyr::distinct() %>%
  dplyr::group_by(Cell, MCP) %>%
  dplyr::filter(abs(MCP_score) > 1) %>%
  dplyr::select(-MCP_score) %>%
  dplyr::distinct()

Metabotype_n_cells = metabotype_meta %>%
  dplyr::filter(Cell %in% MCP_cells_assign$Cell) %>%
  dplyr::group_by(Metabotype) %>%
  dplyr::summarise(n_cells_metabotype = n())

supp_table_21 = MCP_cells_assign %>%
  dplyr::left_join(metabotype_meta) %>%
  dplyr::group_by(MCP, Metabotype) %>%
  dplyr::summarise(n_common_cells = n()) %>%
  dplyr::left_join(Metabotype_n_cells) %>%
  dplyr::mutate(prop_common = n_common_cells / n_cells_metabotype)

write.csv(supp_table_21, file.path(out_dir, "Supp_Table_21_MCP_metabotype_overlap.csv"), row.names = FALSE)


# Tables not referenced in the manuscript (kept for reference) ===================================
# Useful outputs from the same analyses, but not among the 21 numbered supplementary tables.

# Baysor segmentation cell stats
cell_stats = read.csv(file.path(cluster_base_dir, "Results/Baysor_segmentation/Baysor_run_w_segmask_out_dev/segmentation_cell_stats.csv"))

baysor_count_mat = readRDS("../../Data/Figure 1/Baysor_filtered_count_mat.rds")$count_mat
n_genes_per_cell = colSums(!is.na(baysor_count_mat %>% dplyr::select(-gene)))
n_genes_df = data.frame(cell = names(n_genes_per_cell), n_genes = as.numeric(n_genes_per_cell))

extra_table_baysor_cell_stats = cell_stats %>%
  dplyr::select(cell, x, y, n_transcripts) %>%
  dplyr::left_join(n_genes_df, by = "cell")

write.csv(extra_table_baysor_cell_stats, file.path(out_dir, "Extra_Table_Baysor_cell_stats.csv"), row.names = FALSE)


# Region-specific metabolite markers (parent regions)
# Needs a live FindAllMarkers() re-run -- reg_parent_DE isn't persisted anywhere, and
# region_parent isn't a metadata column on the saved combined object.
depth_3_merge = c("amc" = "mfbc", "fa" = "cc", "fxs" = "mfbc")

metabo_RNA_seurat_p12 = metabo_RNA_seurat
a = metabo_RNA_seurat_p12$region
b = lapply(metabo_RNA_seurat_p12$region, get_brain_parent_depth, depth = 4) %>% unlist()
metabo_RNA_seurat_p12@meta.data$region_parent = as.character(b)
metabo_RNA_seurat_p12@meta.data = metabo_RNA_seurat_p12@meta.data %>%
  dplyr::mutate(region_parent = recode(region_parent, !!! depth_3_merge))

cortical_plate_regions = unique(metabo_RNA_seurat_p12@meta.data$region[metabo_RNA_seurat_p12@meta.data$region_parent == "CTXpl"])
d_5_ctxpl = lapply(cortical_plate_regions, get_brain_parent_depth, depth = 5) %>% unlist()
names(cortical_plate_regions) = d_5_ctxpl

a = names(cortical_plate_regions)[match(a, cortical_plate_regions)]
c = a
c[is.na(c)] = b[is.na(c)]
metabo_RNA_seurat_p12@meta.data$region_parent = as.character(c)
metabo_RNA_seurat_p12@meta.data = metabo_RNA_seurat_p12@meta.data %>%
  dplyr::mutate(region_parent = recode(region_parent, !!! depth_3_merge))

reg_parent_DE = FindAllMarkers(object = metabo_RNA_seurat_p12,
                               assay = "Metabo",
                               group.by = "region_parent",
                               only.pos = F,
                               logfc.threshold = 0.25,
                               min.cells.group = 20,
                               min.cells.feature = 20,
                               min.pct = 0.1)

extra_table_region_parent_markers = reg_parent_DE %>%
  dplyr::filter(abs(avg_log2FC) >= 1, p_val_adj < 0.05,
                pct.1 >= 0.5 | pct.2 >= 0.5, cluster != "Uncertain")

write.csv(extra_table_region_parent_markers, file.path(out_dir, "Extra_Table_region_parent_metabolite_markers.csv"), row.names = FALSE)

rm(metabo_RNA_seurat_p12)


# SpatialGlue integration
sel_sg_res = readRDS("../../Data/Figure 6/updated_SG_eval_optim_res.rds")
joint_embed = read.csv("../../Data/Figure 6/updated_Joint_embedding_best_param.csv")

joint_embed_mat = joint_embed
joint_embed_mat$X = sel_sg_res$SG_clusts_metadata$X
colnames(joint_embed_mat) = c("Cell", paste0("SG", c(1:50)))
joint_embed_mat = joint_embed_mat %>% column_to_rownames("Cell")

extra_table_spatialglue_joint_embed = sel_sg_res$SG_clusts_metadata %>%
  dplyr::left_join(joint_embed_mat %>% rownames_to_column("X"), by = "X")

write.csv(extra_table_spatialglue_joint_embed, file.path(out_dir, "Extra_Table_SpatialGlue_joint_embedding.csv"), row.names = FALSE)


# SpatialGlue parameter optimization
all_sg_res = readRDS("../../Data/Figure 6/SG_eval_res_combined.rds")

extra_table_spatialglue_param_optim = all_sg_res$SCIB_metrics %>%
  dplyr::select(-X) %>%
  dplyr::mutate(param_comb = sub(".csv.*", "", param_comb)) %>%
  tidyr::separate(param_comb, into = c("Dim", "K_feat", "metric_feat", "Loss_weights"), sep = "_", remove = F) %>%
  dplyr::mutate(Dim = sub("dim", "", Dim),
                K_feat = sub("k", "", K_feat),
                metric_feat = sub("metric", "", metric_feat),
                Loss_weights = sub("wf", "", Loss_weights)) %>%
  dplyr::mutate(Metric_label = paste0(Metric, ";", Label))

write.csv(extra_table_spatialglue_param_optim, file.path(out_dir, "Extra_Table_SpatialGlue_param_optimization.csv"), row.names = FALSE)

message("Done -- all supplementary tables written to ", out_dir)

# Unnumbered companion export (Compare_DE_DIALOGUE.R)
to_remove_metabo = c("prontosil-H")
MCP_up_markers = MCP_signature_to_df(sig_list = Res$R$sig2, min_celltype = 5) %>%
  dplyr::select(-cell_type) %>%
  dplyr::distinct() %>%
  dplyr::filter(Direction == "up") %>%
  dplyr::left_join(metabo_molec_mapping[, c("ion", "mol_name_one_to_one")], by = c("metabolite" = "ion")) %>%
  dplyr::filter(!mol_name_one_to_one %in% to_remove_metabo) %>%
  dplyr::select(-Direction, -n_CT)

write.csv(MCP_up_markers, file.path(out_dir, "MCP_up_markers_named.csv"), row.names = FALSE)
