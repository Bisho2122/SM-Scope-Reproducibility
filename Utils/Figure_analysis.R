# Analysis/data-wrangling helpers shared across Fig_scripts/Figure_N.R, extracted from
# duplicated or bulky inline blocks so the figure scripts stay focused on plotting.
# Sourced the same way as utils.R/Plotting.R/colors.R.

# Region metadata -----------------------------------------------------------

add_parent_regions = function(seurat_obj){
  depth_3_merge = c("amc" = "mfbc",
                    "fa" = "cc",
                    "fxs" = "mfbc")

  fine_region = seurat_obj$region
  region_parent_depth4 = lapply(seurat_obj$region, get_brain_parent_depth, depth = 4) %>% unlist()

  seurat_obj@meta.data$region_parent = as.character(region_parent_depth4)
  seurat_obj@meta.data = seurat_obj@meta.data %>%
    dplyr::mutate(region_parent = recode(region_parent, !!! depth_3_merge))

  cortical_plate_regions = unique(seurat_obj@meta.data$region[seurat_obj@meta.data$region_parent == "CTXpl"])
  ctxpl_depth5 = lapply(cortical_plate_regions, get_brain_parent_depth, depth = 5) %>% unlist()
  names(cortical_plate_regions) = ctxpl_depth5

  region_parent_final = names(cortical_plate_regions)[match(fine_region, cortical_plate_regions)]
  region_parent_final[is.na(region_parent_final)] = region_parent_depth4[is.na(region_parent_final)]

  seurat_obj@meta.data$region_parent = as.character(region_parent_final)
  seurat_obj@meta.data = seurat_obj@meta.data %>%
    dplyr::mutate(region_parent = recode(region_parent, !!! depth_3_merge))

  return(seurat_obj)
}

# Metabolite naming -----------------------------------------------------------

map_metabolite_names = function(seurat_obj, rename_assay_rownames = TRUE,
                                 mapping_path = "../Data/General/updated_annots_w_1_1_mapping.csv"){
  metabo_molec_mapping = read.csv(mapping_path)
  metabo_molec_mapping = metabo_molec_mapping[metabo_molec_mapping$ion %in% rownames(seurat_obj@assays$Metabo),]

  metabo_molec_mapping$adduct = str_extract(metabo_molec_mapping$ion,
                                            "(\\+|\\-).*")
  metabo_molec_mapping$mol_name_one_to_one = paste0(metabo_molec_mapping$mol_name_one_to_one, metabo_molec_mapping$adduct)

  if (rename_assay_rownames){
    rownames(seurat_obj@assays$Metabo) = metabo_molec_mapping$mol_name_one_to_one[
      metabo_molec_mapping$ion == rownames(seurat_obj@assays$Metabo)]
  }

  return(list(seurat_obj = seurat_obj, mapping = metabo_molec_mapping))
}

# Atlas polygons -----------------------------------------------------------

prep_atlas_polygons = function(atlas_geojson, seurat_obj, brain_ontology){
  affine_mat = rotation(180)
  affine_mat[1,1] = 1
  newly_annot_ft_parts = atlas_geojson$region[atlas_geojson$id == "1009"]

  atlas_geojson %>%
    dplyr::select(-id) %>%
    dplyr::filter(region %in% c(seurat_obj$region, newly_annot_ft_parts)) %>%
    dplyr::mutate(geometry = geometry * affine_mat) %>%
    dplyr::left_join(brain_ontology[,c("acronym", "color_hex_triplet")],
                     by = c("region" = "acronym"))
}

# Dominant category per region -----------------------------------------------------------

# cell_region_category_df must have one row per cell with a region column (region_col)
# and a category column (category_col); both passed as column-name strings so the
# output keeps whatever names the caller already uses downstream.
compute_dominant_category_per_region = function(cell_region_category_df, region_col, category_col){
  total_cells_regions = cell_region_category_df %>%
    dplyr::group_by(across(all_of(region_col))) %>%
    dplyr::summarise(N_total = n())

  region_category_prop = cell_region_category_df %>%
    dplyr::group_by(across(all_of(c(region_col, category_col)))) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::left_join(total_cells_regions, by = region_col) %>%
    dplyr::mutate(prop = n / N_total)

  region_category_prop %>%
    dplyr::group_by(across(all_of(region_col))) %>%
    dplyr::slice_max(prop, n = 1, with_ties = F) %>%
    dplyr::ungroup()
}

# SpatialGlue-leiden dominant-cluster-per-region (Figure_6 only) -- kept separate from
# compute_dominant_category_per_region() above because its denominator is filtered
# (Region != "") while its numerator isn't, an asymmetry specific to this one call site.
compute_dominant_sg_per_region = function(sg_clusts_metadata){
  total_cells_regions = sg_clusts_metadata %>%
    dplyr::select(X, Region) %>%
    dplyr::filter(Region != "") %>%
    dplyr::group_by(Region) %>%
    dplyr::summarise(N_total = n())

  region_SG_prop = sg_clusts_metadata %>%
    dplyr::select(X,Region,SG_leiden) %>%
    dplyr::group_by(Region, SG_leiden) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::left_join(total_cells_regions) %>%
    dplyr::mutate(prop = n / N_total) %>%
    dplyr::mutate(SG_leiden = paste0("SG_", SG_leiden))

  region_SG_prop %>%
    dplyr::group_by(Region) %>%
    dplyr::slice_max(prop, n = 1, with_ties = F) %>%
    dplyr::ungroup()
}

# DE marker filtering -----------------------------------------------------------

filter_de_markers = function(de_df, lfc_cutoff = 1, lfc_strict = TRUE, padj_cutoff = 0.05,
                              pct_cutoff = 0.5, directional = FALSE,
                              cluster_exclude = NULL, gene_exclude = NULL){
  result = de_df

  if (lfc_strict){
    result = result %>% dplyr::filter(abs(avg_log2FC) > lfc_cutoff)
  } else {
    result = result %>% dplyr::filter(abs(avg_log2FC) >= lfc_cutoff)
  }

  result = result %>% dplyr::filter(p_val_adj < padj_cutoff)

  if (directional){
    result = result %>%
      dplyr::filter((avg_log2FC > 0 & pct.1 >= pct_cutoff) | (avg_log2FC < 0 & pct.2 >= pct_cutoff))
  } else {
    result = result %>%
      dplyr::filter(pct.1 >= pct_cutoff | pct.2 >= pct_cutoff)
  }

  if (!is.null(cluster_exclude)){
    result = result %>% dplyr::filter(!cluster %in% cluster_exclude)
  }
  if (!is.null(gene_exclude)){
    result = result %>% dplyr::filter(!gene %in% gene_exclude)
  }

  return(result)
}

# DE dotplot (aggregate + plot_DE_clusterd_dotplot + common styling) --------------------

de_dotplot = function(seurat_obj, marker_df, assay, group_by, rotate_x = FALSE){
  agg = AverageExpression(seurat_obj, assays = assay, group.by = group_by, slot = "data")[[1]]

  p = plot_DE_clusterd_dotplot(marker_df = marker_df, agg_scaled_mat = agg) +
    scale_color_gradient2(
      low = scale_cols$div$low,
      mid = scale_cols$div$mid,
      high = scale_cols$div$high,
      midpoint = 0
    ) +
    labs(size = "Percentage\nDetection", color = "Z-Score") +
    theme(text = element_text(size = 16))

  if (rotate_x){
    p = p + rotate_x_text(angle = 90)
  }

  return(p)
}

# hclust-based row/column reordering for a long-format df's factor levels --------------------

reorder_by_hclust = function(df, row_var, col_var, value_var, fill = 0){
  mat = df %>%
    dplyr::select(all_of(c(row_var, col_var, value_var))) %>%
    dplyr::distinct() %>%
    pivot_wider(names_from = all_of(col_var), values_from = all_of(value_var), values_fill = fill) %>%
    column_to_rownames(row_var) %>%
    as.matrix()

  row_order = hclust(dist(mat))$order
  col_order = hclust(dist(t(mat)))$order

  df[[row_var]] = factor(df[[row_var]], levels = rownames(mat)[row_order])
  df[[col_var]] = factor(df[[col_var]], levels = colnames(mat)[col_order])

  return(df)
}

# Adjusted Mutual Information vs. region -----------------------------------------------------------

compute_ami_by_region = function(seurat_obj, include_sg = FALSE){
  AMI_data = seurat_obj@meta.data %>%
    rownames_to_column("Cell") %>%
    dplyr::select(Cell, region, Metabotype, Cell_type_fine, SG_leiden) %>%
    dplyr::distinct() %>%
    dplyr::filter(region != "Uncertain") %>%
    dplyr::mutate(CT_coarse = sub("_.*","",Cell_type_fine))

  AMI_CT_fine = aricode::AMI(AMI_data$region, AMI_data$Cell_type_fine)
  AMI_metabotype = aricode::AMI(AMI_data$region, AMI_data$Metabotype)
  AMI_CT_coarse = aricode::AMI(AMI_data$region, AMI_data$CT_coarse)

  clustering_labels = c("Cell type\ncoarse", "Cell type\nfine", "Metabotype")
  ami_values = c(AMI_CT_coarse, AMI_CT_fine, AMI_metabotype)
  result = list(AMI_data = AMI_data, AMI_CT_fine = AMI_CT_fine,
                AMI_metabotype = AMI_metabotype, AMI_CT_coarse = AMI_CT_coarse)

  if (include_sg){
    AMI_SG = aricode::AMI(AMI_data$region, AMI_data$SG_leiden)
    clustering_labels = c(clustering_labels, "Spatial\nGlue")
    ami_values = c(ami_values, AMI_SG)
    result$AMI_SG = AMI_SG
  }

  result$plot_data = data.frame(Clustering = clustering_labels, AMI = ami_values)
  return(result)
}

plot_ami_barplot = function(ami_plot_data){
  ami_plot_data %>%
    ggplot(aes(x = Clustering, y = AMI)) +
    geom_bar(stat = "identity",
             fill = my_cols$highlight_colors[1],
             width = 0.7) +
    theme_pubr() +
    xlab("") +
    ylab("Adjusted Mutual Information") +
    theme(text = element_text(size = 16)) +
    scale_x_discrete(expand = c(0.15,0.15))
}

# Group composition percentages (for stacked-bar plots) -----------------------------------------------------------

compute_group_percent = function(df, count_vars, percent_group_var){
  df %>%
    dplyr::count(across(all_of(count_vars))) %>%
    dplyr::group_by(across(all_of(percent_group_var))) %>%
    dplyr::mutate(Percent = n / sum(n) * 100)
}

# DIALOGUE MCP cell embedding -----------------------------------------------------------

get_mcp_cell_embed = function(dialogue_res){
  dialogue_res$full_RA_res$scores %>%
    dplyr::bind_rows() %>%
    rownames_to_column("Cell")
}

# Figure_4-specific: gene-metabolite / KEGG pathway heatmap matrices --------------------------

build_gene_metabo_pathway_mats = function(niches_corr, sf_name_unique, all_dists, kegg_mouse_pathways){
  niches_overlap_rxn = niches_corr %>%
    dplyr::filter(abs(coef) >= 0.5,
                  n_cells_cor >= 30) %>%
    dplyr::mutate(feat_A = sub(".*rna_", "", feat_A),
                  feat_B = sub(".*metabo_", "", feat_B),
                  feat_B = sub("[+-].*", "", feat_B)) %>%
    dplyr::left_join(all_dists,by = c("feat_A" = "gene",
                                      "feat_B" = "metabo")) %>%
    dplyr::filter(!is.na(pathway_name))

  best_dist_per_pair = niches_overlap_rxn %>%
    dplyr::select(Niche, feat_A, feat_B,coef,
                  pathway_source, pathway_name,dist) %>%
    dplyr::distinct() %>%
    dplyr::group_by(Niche, feat_A, feat_B, pathway_name) %>%
    dplyr::slice_min(dist, n = 1, with_ties = F)

  cor_info_mat = best_dist_per_pair %>%
    ungroup() %>%
    dplyr::select(Niche, feat_A, feat_B, coef) %>%
    dplyr::left_join(sf_name_unique, by = c("feat_B" = "ion")) %>%
    dplyr::select(-feat_B) %>%
    dplyr::rename("feat_B" = "mol_name_one_to_one") %>%
    dplyr::distinct() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(gene_metabo = paste0(feat_A, " - ", feat_B)) %>%
    dplyr::select(gene_metabo, Niche, coef) %>%
    pivot_wider(names_from = Niche, values_from = coef, values_fill = 0) %>%
    column_to_rownames("gene_metabo") %>%
    as.matrix()

  pathway_info_mat = best_dist_per_pair %>%
    ungroup() %>%
    dplyr::select(feat_A, feat_B, pathway_source, pathway_name, dist) %>%
    dplyr::left_join(sf_name_unique, by = c("feat_B" = "ion")) %>%
    dplyr::select(-feat_B) %>%
    dplyr::rename("feat_B" = "mol_name_one_to_one") %>%
    dplyr::distinct() %>%
    dplyr::mutate(pathway_name = gsub("rn", "mmu", pathway_name)) %>%
    dplyr::left_join(kegg_mouse_pathways, by = c("pathway_name" = "kegg_pathway_id")) %>%
    dplyr::mutate(kegg_pathway_name = sub(" - Mus musculus.*", "", kegg_pathway_name)) %>%
    dplyr::mutate(kegg_pathway_name = ifelse(is.na(kegg_pathway_name),
                                             pathway_name,
                                             kegg_pathway_name)) %>%
    dplyr::select(-pathway_name) %>%
    dplyr::rename("pathway_name" = "kegg_pathway_name") %>%
    dplyr::distinct() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(gene_metabo = paste0(feat_A, " - ", feat_B)) %>%
    dplyr::select(gene_metabo, pathway_name) %>%
    dplyr::mutate(in_path = 1) %>%
    pivot_wider(names_from = pathway_name, values_from = in_path, values_fill = 0) %>%
    column_to_rownames("gene_metabo") %>%
    as.matrix()

  pathway_info_label_mat = best_dist_per_pair %>%
    ungroup() %>%
    dplyr::select(feat_A, feat_B, pathway_source, pathway_name, dist) %>%
    dplyr::left_join(sf_name_unique, by = c("feat_B" = "ion")) %>%
    dplyr::select(-feat_B) %>%
    dplyr::rename("feat_B" = "mol_name_one_to_one") %>%
    dplyr::distinct() %>%
    dplyr::mutate(pathway_name = gsub("rn", "mmu", pathway_name)) %>%
    dplyr::left_join(kegg_mouse_pathways, by = c("pathway_name" = "kegg_pathway_id")) %>%
    dplyr::mutate(kegg_pathway_name = sub(" - Mus musculus.*", "", kegg_pathway_name)) %>%
    dplyr::mutate(kegg_pathway_name = ifelse(is.na(kegg_pathway_name),
                                             pathway_name,
                                             kegg_pathway_name)) %>%
    dplyr::select(-pathway_name) %>%
    dplyr::rename("pathway_name" = "kegg_pathway_name") %>%
    dplyr::distinct() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(gene_metabo = paste0(feat_A, " - ", feat_B)) %>%
    dplyr::select(gene_metabo, pathway_name, dist) %>%
    pivot_wider(names_from = pathway_name, values_from = dist, values_fill = -1) %>%
    column_to_rownames("gene_metabo") %>%
    as.matrix()

  pathway_info_mat = pathway_info_mat[rownames(cor_info_mat),]
  pathway_info_label_mat = pathway_info_label_mat[rownames(cor_info_mat),
                                                  colnames(pathway_info_mat)]
  pathway_info_label_mat[pathway_info_label_mat == -1] = ""

  return(list(cor_info_mat = cor_info_mat,
              pathway_info_mat = pathway_info_mat,
              pathway_info_label_mat = pathway_info_label_mat,
              niches_of_interest = unique(best_dist_per_pair$Niche)))
}

# Figure_4-specific: with/without-intra MISTY RMSE comparison panels --------------------------

misty_rmse_comparison_plots = function(with_intra_dir, without_intra_dir){
  with_intra_misty_res = mistyR::collect_results(with_intra_dir)
  without_intra_misty_res = mistyR::collect_results(without_intra_dir)

  with_intra_RMSE = with_intra_misty_res$improvements %>%
    dplyr::filter(measure %in% c("gain.RMSE", "intra.RMSE",
                                 "multi.RMSE", "p.RMSE")) %>%
    dplyr::mutate(measure = paste0("intra", "_", measure)) %>%
    spread(key = "measure", value = "value")

  without_intra_RMSE = without_intra_misty_res$improvements %>%
    dplyr::filter(measure %in% c("gain.RMSE", "intra.RMSE",
                                 "multi.RMSE", "p.RMSE")) %>%
    dplyr::mutate(measure = paste0("no.intra", "_", measure)) %>%
    spread(key = "measure", value = "value")

  rmse_comparison = cbind.data.frame(with_intra_RMSE, without_intra_RMSE) %>%
    dplyr::select(-target, -sample)

  p1 = rmse_comparison %>%
    ggplot(aes(x = intra_multi.RMSE, y = no.intra_multi.RMSE)) +
    geom_point(size = 2) +
    geom_abline() +
    theme_pubr()

  p2 = rmse_comparison %>%
    ggplot(aes(x = intra_intra.RMSE, y = no.intra_intra.RMSE)) +
    geom_point(size = 2) +
    geom_abline() +
    theme_pubr()

  p3 = rmse_comparison %>%
    ggplot(aes(x = intra_gain.RMSE, y = no.intra_gain.RMSE,
               colour = intra_gain.RMSE > 1)) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 1) +
    scale_color_manual(values = assign_colors(c("FALSE", "TRUE"),
                                              my_cols$binary)) +
    theme_pubr()

  # trace edited plot_view_contributions to return a ggplot object instead of printing
  p4 = mistyR::plot_view_contributions(with_intra_misty_res,trim.measure = "gain.RMSE",
                                                        trim = 1) +
    geom_hline(yintercept = 0.75, linetype = "dashed") +
    scale_fill_manual(values = assign_colors(c("gene", "intra"),
                                             my_cols$highlight_colors))

  ggarrange(p2, p1, p3, p4, labels = c("A", "B", "C", "D"))
}

# Per-figure manuscript-text number exports ------------------------------------
# One function per figure, called by Fig_scripts/Generate_Figure_Numbers.R. Each takes the
# already-loaded input objects for that figure and reconstructs whatever intermediate values
# its numbers need, reusing the shared helpers above where the same computation is also used
# for a plot.

compute_figure2_numbers = function(metabo_RNA_seurat, DE_res){
  metabo_RNA_seurat@meta.data$celltype_coarse = sub("_.*", "", metabo_RNA_seurat@meta.data$Cell_type_fine)

  n_metabo_types = length(unique(metabo_RNA_seurat$Metabotype))
  n_celltypes = length(unique(metabo_RNA_seurat$Cell_type_fine))
  n_celltypes_coarse = length(unique(metabo_RNA_seurat$celltype_coarse))
  n_rna_leiden = length(unique(metabo_RNA_seurat$seurat_clusters))

  n_de_metabo_celltype_markers = nrow(filter_de_markers(DE_res$Metabo$Celltypes))
  n_de_rna_metabotype_markers = nrow(filter_de_markers(DE_res$RNA$Metabotype))
  n_de_metabo_metabotype_markers = nrow(filter_de_markers(DE_res$Metabo$Metabotype))

  data.frame(
    metric = c("n_rna_leiden_clusters", "n_celltypes_fine", "n_celltypes_coarse",
               "n_metabotypes", "n_metabolites", "n_genes",
               "n_de_metabo_celltype_markers", "n_de_rna_metabotype_markers",
               "n_de_metabo_metabotype_markers"),
    value = c(n_rna_leiden, n_celltypes, n_celltypes_coarse,
              n_metabo_types, nrow(metabo_RNA_seurat@assays$Metabo), nrow(metabo_RNA_seurat@assays$RNA),
              n_de_metabo_celltype_markers, n_de_rna_metabotype_markers,
              n_de_metabo_metabotype_markers)
  )
}

compute_figure3_numbers = function(metabo_RNA_seurat, SNR_fisher_res){
  metabo_RNA_seurat = add_parent_regions(metabo_RNA_seurat)

  n_regions_with_misty_important = nrow(SNR_fisher_res)
  regions_with_high_SNR_enrich = sum(SNR_fisher_res$fisher_pval_log > -log10(0.05))

  data.frame(
    metric = c("n_regions_fine", "n_regions_parent",
               "n_regions_with_misty_important", "regions_with_high_SNR_enrich"),
    value = c(length(unique(metabo_RNA_seurat$region)),
              length(unique(metabo_RNA_seurat$region_parent)),
              n_regions_with_misty_important, regions_with_high_SNR_enrich)
  )
}

# Returns a list: $numbers (figure_4_numbers.csv) and $n_cells_in_VAL_example
# (figure_4_VAL_example_n_cells_per_celltype.csv).
compute_figure4_numbers = function(all_corr_res, all_dists, brain_ontology, metabo_molec_mapping){
  Thalamas_region_only = brain_ontology$acronym[
    str_detect(brain_ontology$structure_id_path, "549/")
    ]

  cor_CT_regions = all_corr_res$Cell_type_and_Region %>%
    tidyr::separate(CT_Region, into = c("Cell_type", "Region"), sep = "[.]", remove = F) %>%
    dplyr::mutate(feat_A = sub(".*rna_", "", feat_A),
                  feat_B = sub(".*metabo_", "", feat_B)) %>%
    dplyr::filter(abs(coef) > 0.5,
                  n_cells_cor >= 30,
                  Region %in% Thalamas_region_only) %>%
    dplyr::filter(feat_B %in% metabo_molec_mapping$mol_name_one_to_one)

  n_gene_metabo_pairs = all_corr_res$All_cells %>%
    dplyr::select(feat_A, feat_B) %>%
    dplyr::distinct() %>%
    nrow()
  n_pairs_abs_0_5 = all_corr_res$All_cells %>%
    dplyr::filter(abs(coef) >= 0.5) %>%
    dplyr::select(feat_A, feat_B) %>%
    dplyr::distinct() %>%
    nrow()
  n_pairs_pathway_mapped = nrow(all_dists)
  n_cells_in_VAL_example = cor_CT_regions %>%
    dplyr::filter(Region == "VAL") %>%
    dplyr::select(Cell_type, n_cells_cor) %>%
    dplyr::distinct()

  list(
    numbers = data.frame(
      metric = c("n_gene_metabo_pairs", "n_pairs_abs_0_5", "n_pairs_pathway_mapped"),
      value = c(n_gene_metabo_pairs, n_pairs_abs_0_5, n_pairs_pathway_mapped)
    ),
    n_cells_in_VAL_example = n_cells_in_VAL_example
  )
}

compute_figure5_numbers = function(metabo_RNA_seurat, dialogue_res, all_DE_res, spatial_niches){
  n_niches_region = spatial_niches %>%
    dplyr::select(Region, Niche) %>%
    dplyr::distinct() %>%
    dplyr::count(Region)

  dialogue_sig_min5 = MCP_signature_to_df(sig_list = dialogue_res$R$sig2, min_celltype = 5) %>%
    dplyr::select(-cell_type) %>%
    dplyr::distinct()
  dialogue_sig_min4 = MCP_signature_to_df(sig_list = dialogue_res$R$sig2, min_celltype = 4) %>%
    dplyr::select(-cell_type) %>%
    dplyr::distinct()

  mcp2_metabos_min5 = unique(dialogue_sig_min5$metabolite[dialogue_sig_min5$MCP == "MCP2"])
  other_mcp_metabos_min5 = unique(dialogue_sig_min5$metabolite[dialogue_sig_min5$MCP != "MCP2"])
  mcp2_metabos_min4 = unique(dialogue_sig_min4$metabolite[dialogue_sig_min4$MCP == "MCP2"])

  dialogue_res$R$param$averaging.function = colMeans2
  ct_corr_matrix = DIALOGUE:::DIALOGUE.plot.av(dialogue_res$R, MCPs = 1:10)
  mcp2_corr_ct = cor(ct_corr_matrix[[2]])
  mcp2_corr_lower_tri = mcp2_corr_ct[lower.tri(mcp2_corr_ct)]

  Met1_markers = all_DE_res$Metabo$Metabotype %>%
    dplyr::filter(avg_log2FC >= 1, p_val_adj < 0.05, cluster == "Met.1") %>%
    dplyr::pull(gene)

  MCP_cell_embed = get_mcp_cell_embed(dialogue_res)
  metabotype_meta = metabo_RNA_seurat@meta.data %>%
    rownames_to_column("Cell") %>%
    dplyr::select(Cell, Metabotype)
  conf_mat_MCP_metabotype = MCP_cell_embed %>%
    dplyr::left_join(metabotype_meta)
  conf_mat_MCP_metabotype = conf_mat_MCP_metabotype[,c("Cell",
                                                       paste0("MCP", 1:10),
                                                       "Metabotype")]
  conf_mat_MCP_metabotype = conf_mat_MCP_metabotype %>%
    gather(key = "MCP", value = "MCP_score", -Cell, -Metabotype)

  MCP_cells_assign = conf_mat_MCP_metabotype %>%
    dplyr::select(Cell, MCP, MCP_score) %>%
    dplyr::distinct() %>%
    dplyr::group_by(Cell,MCP) %>%
    dplyr::filter(abs(MCP_score) > 1) %>%
    dplyr::select(-MCP_score) %>%
    dplyr::distinct()

  Metabotype_n_cells = metabotype_meta %>%
    dplyr::filter(Cell %in% MCP_cells_assign$Cell) %>%
    dplyr::group_by(Metabotype) %>%
    dplyr::summarise(n_cells_metabotype = n())

  conf_mat_MCP_metabotype = MCP_cells_assign %>%
    dplyr::left_join(metabotype_meta) %>%
    dplyr::group_by(MCP, Metabotype) %>%
    dplyr::summarise(n_common_cells = n()) %>%
    dplyr::left_join(Metabotype_n_cells) %>%
    dplyr::mutate(prop_common = n_common_cells / n_cells_metabotype)

  mcp2_met1_pct_overlap = conf_mat_MCP_metabotype %>%
    dplyr::filter(MCP == "MCP2", Metabotype == "Met.1") %>%
    dplyr::pull(prop_common)

  data.frame(
    metric = c("n_spatial_niches", "n_regions_with_niches", "median_niches_per_region",
               "n_mcps", "mcp2_n_metabolites_min_celltype_5", "mcp2_n_metabolites_min_celltype_4",
               "mcp2_n_unique_metabolites_min_celltype_5", "mcp2_min_celltype_corr",
               "met1_n_de_markers", "mcp2_met1_n_overlap", "mcp2_met1_pct_overlap"),
    value = c(length(unique(spatial_niches$Niche)),
              length(unique(spatial_niches$Region)),
              median(n_niches_region$n),
              10,
              length(mcp2_metabos_min5),
              length(mcp2_metabos_min4),
              length(setdiff(mcp2_metabos_min5, other_mcp_metabos_min5)),
              min(mcp2_corr_lower_tri),
              length(Met1_markers),
              length(intersect(mcp2_metabos_min5, Met1_markers)),
              mcp2_met1_pct_overlap)
  )
}

# Returns a list: $numbers (figure_6_numbers.csv) and $dominant_SG_per_region
# (figure_6_dominant_SG_per_region.csv). `sel_sg_res` must already have the Cell_type_fine/coarse
# relabeling applied (see Fig_scripts/Figure_6.R's "Subset data for selected param comb" section).
compute_figure6_numbers = function(metabo_RNA_seurat, sel_sg_res){
  metabo_RNA_seurat = add_parent_regions(metabo_RNA_seurat)

  cell_region_parent = metabo_RNA_seurat@meta.data %>%
    rownames_to_column("Cell") %>%
    dplyr::select(Cell, region_parent) %>%
    dplyr::distinct()

  sg_leiden_filtered = sel_sg_res$SG_clusts_metadata %>%
    dplyr::mutate(SG_leiden = paste0("SG_", SG_leiden)) %>%
    dplyr::rename(Cell = X) %>%
    dplyr::left_join(cell_region_parent, by = "Cell") %>%
    dplyr::filter(Cell_type_coarse != "", region_parent != "Uncertain")

  ami_res = compute_ami_by_region(metabo_RNA_seurat, include_sg = TRUE)
  region_SG_prop = compute_dominant_sg_per_region(sel_sg_res$SG_clusts_metadata)

  numbers = data.frame(
    metric = c("n_sg_leiden_clusters", "ami_sg_vs_region", "ami_metabotype_vs_region",
               "ami_celltype_fine_vs_region", "ami_celltype_coarse_vs_region",
               "mean_alpha_metabo", "mean_alpha_rna",
               "n_regions_with_dominant_sg_over_65pct"),
    value = c(length(unique(sg_leiden_filtered$SG_leiden)),
              ami_res$AMI_SG, ami_res$AMI_metabotype, ami_res$AMI_CT_fine, ami_res$AMI_CT_coarse,
              mean(sel_sg_res$SG_clusts_metadata$alpha_metabo, na.rm = TRUE),
              mean(sel_sg_res$SG_clusts_metadata$alpha_RNA, na.rm = TRUE),
              sum(region_SG_prop$prop > 0.65))
  )

  list(numbers = numbers,
       dominant_SG_per_region = region_SG_prop %>% dplyr::select(Region, SG_leiden, prop))
}
