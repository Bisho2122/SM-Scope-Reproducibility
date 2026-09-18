library(fplot)

# General -----------------------------------------------------------------
setFplot_page(page = "us", margins = 1)
export_plot = function(plot, path, ...){
  export_graph_start(file = path, ...)
  print(plot)
  export_graph_end()
}

clean_labels <- function(p) {
  # --- helpers ---
  clean_text <- function(txt) {
    if (is.null(txt)) return(NULL)
    txt <- str_replace_all(txt, "_", " ")
    if (str_count(txt, "\\s+") >= 1) txt <- str_to_sentence(txt)
    txt
  }
  clean_vec <- function(x) vapply(x, clean_text, character(1))

  # --- 1) Clean plot/legend titles ---
  labs_list <- lapply(p$labels, clean_text)
  p <- p + do.call(labs, labs_list)

  # --- 2) Inspect built plot to detect discrete vs continuous ---
  built <- ggplot_build(p)
  panel_params <- built$layout$panel_params[[1]]

  x_is_discrete <- inherits(panel_params$x$scale, "ScaleDiscrete")
  y_is_discrete <- inherits(panel_params$y$scale, "ScaleDiscrete")

  if (x_is_discrete) p <- p + scale_x_discrete(labels = clean_vec)
  if (y_is_discrete) p <- p + scale_y_discrete(labels = clean_vec)

  # --- 3) Legend entries: handle only discrete aesthetics ---
  for (aes in c("fill", "colour", "color", "shape", "linetype")) {
    scale <- built$plot$scales$get_scales(aes)
    if (!is.null(scale) && inherits(scale, "ScaleDiscrete")) {
      if (aes %in% c("colour", "color")) {
        p <- p + scale_color_discrete(labels = clean_vec)
      } else if (aes == "fill") {
        p <- p + scale_fill_discrete(labels = clean_vec)
      } else if (aes == "shape") {
        p <- p + scale_shape_discrete(labels = clean_vec)
      } else if (aes == "linetype") {
        p <- p + scale_linetype_discrete(labels = clean_vec)
      }
    }
  }

  # --- 4) Facet strip labels ---
  if (!inherits(p$facet, "FacetNull")) {
    old_labeller <- p$facet$params$labeller
    if (is.null(old_labeller)) old_labeller <- label_value
    p$facet$params$labeller <- function(labels) {
      labs_df <- old_labeller(labels)
      as.data.frame(lapply(labs_df, clean_vec))
    }
  }

  p
}
# helper to test if a scale is discrete
is.discrete <- function(scale) {
  if (is.null(scale)) return(FALSE)
  inherits(scale, "ScaleDiscrete")
}

add_small_umap_axes <- function(
    p,
    reduction       = "umap",
    corner_offset   = 0.03,
    axis_fraction   = 0.2,
    lineTextcol     = "black",
    cornerTextSize  = 3
) {
  gb <- ggplot_build(p)
  d  <- gb$data[[1]]

  xmin <- min(d$x, na.rm = TRUE)
  xmax <- max(d$x, na.rm = TRUE)
  ymin <- min(d$y, na.rm = TRUE)
  ymax <- max(d$y, na.rm = TRUE)

  x_range <- xmax - xmin
  y_range <- ymax - ymin

  # Origin inside the corner
  x0 <- xmin - corner_offset * x_range
  y0 <- ymin - corner_offset * y_range

  # Axis endpoints
  x1 <- x0 + axis_fraction * x_range
  y1 <- y0 + axis_fraction * y_range

  # Midpoints
  x_mid <- (x0 + x1) / 2
  y_mid <- (y0 + y1) / 2

  # Label text
  if (startsWith(reduction, "umap")) {
    axs_label <- c("UMAP2", "UMAP1")
  } else if (startsWith(reduction, "tsne")) {
    axs_label <- c("t-SNE2", "t-SNE1")
  } else {
    stop("reduction must start with 'umap' or 'tsne'")
  }

  # Label offsets (small shifts)
  x_offset <- 0.03 * y_range  # vertical shift for x-axis label
  y_offset <- 0.03 * x_range  # horizontal shift for y-axis label

  # Axis segments (both start at same origin)
  axes_df <- data.frame(
    x    = c(x0, x0),
    xend = c(x0, x1),
    y    = c(y0, y0),
    yend = c(y1, y0)
  )

  # Correctly shifted label positions
  label_df <- data.frame(
    lab   = axs_label,
    angle = c(90, 0),  # y-axis vertical, x-axis horizontal
    x     = c(x0 - y_offset, x_mid),
    y     = c(y_mid,        y0 - x_offset)
  )

  p +
    geom_segment(
      data        = axes_df,
      inherit.aes = FALSE,
      aes(x = x, xend = xend, y = y, yend = yend),
      arrow = arrow(length = unit(0.1, "inches"), type = "closed"),
      color = lineTextcol
    ) +
    geom_text(
      data        = label_df,
      inherit.aes = FALSE,
      aes(x = x, y = y, label = lab, angle = angle),
      color    = lineTextcol,
      size     = cornerTextSize,
      fontface = "italic"
    )
}

# Regions - Section 3 ------------------------------------------------------
Plot_confusion_hm = function(seurat_obj, clust_1_column,
                             clust_2_column, reference_column,
                             pct_thresh = 1){
  confusion_metadata = seurat_obj@meta.data %>%
    rownames_to_column("Cell") %>%
    dplyr::select(Cell, {{clust_1_column}}, {{clust_2_column}})

  conf_mat <- table(confusion_metadata[[clust_1_column]],
                    confusion_metadata[[clust_2_column]])

  if(reference_column == clust_1_column){
    conf_mat_norm <- prop.table(conf_mat, margin = 1) * 100
  }
  else if (reference_column == clust_2_column){
    conf_mat_norm <- prop.table(conf_mat, margin = 2) * 100
  }
  else{
    stop("Reference column has to match with one of cluster columns")
  }

  # # Convert to long format
  conf_df <- as.data.frame(conf_mat_norm) %>%
    rename(!!clust_1_column := Var1,
           !!clust_2_column := Var2,
           Percentage = Freq)

  #Plot heatmap
  p = ggplot(conf_df, aes_string(x = clust_1_column,
                                 y = clust_2_column, fill = "Percentage")) +
    geom_tile() +
    scale_fill_gradient2(low = "white", high = "red") +
    theme_minimal() +
    rotate_x_text(angle = 45) +
    labs(title = paste0(clust_1_column, " vs ", clust_2_column),
         x = clust_1_column, y = clust_2_column, fill = "Percentage") +
    # theme(text = element_text(size = 12)) +
    geom_text(aes(label = ifelse(Percentage < pct_thresh, "", sprintf("%.1f%%", Percentage))),
              color = "black", size = 3)
  return(p)
}

plot_DE_hm = function(marker_df,
                      LFC_cutoff = 1,
                      padj_cutoff = 0.05,
                      min_pct = 0.5,
                      font_size = 12,
                      top_x = 3){
  filt_df = marker_df %>%
    dplyr::filter(avg_log2FC >= LFC_cutoff,
                  p_val_adj < padj_cutoff,
                  pct.1 >= min_pct) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(avg_log2FC,n = top_x)

  hm_mat = filt_df %>%
    dplyr::select(gene, cluster, avg_log2FC) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(values_from = avg_log2FC, names_from = cluster) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()

  hm_mat[is.na(hm_mat)] = 0
  cols = pals::brewer.blues(n = 50)
  cols[1] = "#FFFFFF"
  pheatmap::pheatmap(hm_mat,color = cols,
           cluster_rows = T,fontsize = font_size)
}

plot_DE_clusterd_dotplot = function(marker_df,
                                    agg_scaled_mat){
  agg_long = agg_scaled_mat %>%
    as.matrix() %>%
    t() %>% scale() %>% t() %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cluster", values_to = "Z_score")

  agg_long$cluster = gsub("-", "_", agg_long$cluster)
  marker_df$cluster = gsub("-", "_", marker_df$cluster)

  marker_df_scored = marker_df %>%
    dplyr::left_join(agg_long, by = c("cluster","gene"))

  mat <- marker_df_scored %>%
    dplyr::select(gene, cluster, Z_score) %>%
    dplyr::distinct() %>%
    pivot_wider(names_from = cluster, values_from = Z_score, values_fill = 0) %>%
    column_to_rownames("gene")

  gene_dist <- dist(mat)
  gene_clust <- hclust(gene_dist)
  gene_order <- gene_clust$labels[gene_clust$order]

  cluster_dist <- dist(t(mat))
  cluster_clust <- hclust(cluster_dist)
  cluster_order <- cluster_clust$labels[cluster_clust$order]

  marker_df_scored$gene <- factor(marker_df_scored$gene, levels = gene_order)
  marker_df_scored$cluster <- factor(marker_df_scored$cluster, levels = cluster_order)

  p_dot = ggplot(marker_df_scored, aes(x = cluster, y = gene)) +
    geom_point(aes(size = pct.1, color = Z_score)) +
    scale_color_gradient2(low = "darkblue", high = "darkred",mid = "white") +
    scale_size(range = c(1, 4)) +
    theme_bw() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(p_dot)

}


# Metabo-gene corr - Section 4 ---------------------------------------------

plot_spatial_cor_res = function(seurat_obj,Region_name,
                                cor_input_data,
                                rna_feat,
                                metabo_feat,
                                coef,
                                cell_type = NULL,
                                convert_sf_to_name = F){

  seurat_obj@graphs = list()
  reg_cells = cor_input_data$Region[[Region_name]] %>% rownames()
  reg_seurat = SeuratObject:::subset.Seurat(seurat_obj,
                                            cells = reg_cells)

  p1 = SpatialDimPlot(seurat_obj,combine = F,
                      cells.highlight = list("Cells in Region" = reg_cells),
                      cols.highlight = c("red", "lightgrey"),
                      alpha = 1,
                      stroke = NA)[[1]] +
    ggtitle(Region_name,
            subtitle = paste0("N cells = ", length(reg_cells))) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

  if(!is.null(cell_type)){
    ct_reg_name = paste0(cell_type,
                         ".",
                         Region_name)
    reg_ct_cells = cor_input_data$Cell_type_and_Region[[ct_reg_name]] %>%
      rownames()

    reg_ct_seurat = SeuratObject:::subset.Seurat(seurat_obj,
                                                 cells = reg_ct_cells)

    reg_corr_cells = get_cor_cells_input(impute_ind_mat = cor_input_data$Cell_type_and_Region[[ct_reg_name]],
                                         rna_feat = rna_feat,
                                         metabo_feat = metabo_feat)

    cell_highlight_list_ct = list(reg_ct_cells)
    names(cell_highlight_list_ct)[1] = cell_type

    p2 = SpatialDimPlot(reg_seurat,combine = F,
                        pt.size.factor = 5,
                        cells.highlight = cell_highlight_list_ct,
                        cols.highlight = c("red", "lightgrey"),
                        stroke = NA,
                        alpha = 1)[[1]] +
      ggtitle(Region_name,
              subtitle = paste0("N cells highlighted = ", length(reg_ct_cells))) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
    p3 = SpatialDimPlot(reg_ct_seurat,combine = F,
                        pt.size.factor = 5,
                        cells.highlight = list("Cells for correlation" = reg_corr_cells),
                        cols.highlight = c("red", "lightgrey"),
                        stroke = NA,
                        alpha = 1)[[1]] +
      ggtitle(ct_reg_name,
              subtitle = paste0("N cells highlighted = ", length(reg_corr_cells))) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  }
  else{
    reg_corr_cells = get_cor_cells_input(impute_ind_mat = cor_input_data$Region[[Region_name]],
                                         rna_feat = rna_feat,
                                         metabo_feat = metabo_feat)

    p2 = SpatialDimPlot(reg_seurat,combine = F,
                        pt.size.factor = 5,
                        cells.highlight = list("Cells for correlation" = reg_corr_cells),
                        cols.highlight = c("red", "lightgrey"),
                        stroke = NA,
                        alpha = 1)[[1]] +
      ggtitle(Region_name,
              subtitle = paste0("N cells highlighted = ", length(reg_corr_cells))) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  }

  reg_corr_seurat = SeuratObject:::subset.Seurat(seurat_obj,
                                                 cells = reg_corr_cells)
  DefaultAssay(reg_corr_seurat) = "RNA"
  feat_A_name = sub(".*rna_","",rna_feat)
  feat_A_p = SpatialPlot(object = reg_corr_seurat,
                         features = feat_A_name,
                         combine = F,
                         pt.size.factor = 8)

  DefaultAssay(reg_corr_seurat) = "Metabo"
  if(convert_sf_to_name){
    feat_B_name = sf_to_name(sub(".*metabo_","",metabo_feat))
  }
  else{
    feat_B_name = sub(".*metabo_","",metabo_feat)
  }
  feat_B_p = SpatialPlot(object = reg_corr_seurat,
                         features = feat_B_name,
                         combine = F,
                         pt.size.factor = 8,slot = "data")

  feat_long = feat_A_p[[1]]$data[,c("x","y","cells",feat_A_name)] %>%
    dplyr::left_join(feat_B_p[[1]]$data[,c("x","y","cells",feat_B_name)],
                     by = c("x","y","cells")) %>%
    gather(key = "Feature", value = "Log_Expression", -x, -y, -cells)

  scatter_data = feat_A_p[[1]]$data[,c("x","y","cells",feat_A_name)] %>%
    dplyr::left_join(feat_B_p[[1]]$data[,c("x","y","cells",feat_B_name)],
                     by = c("x","y","cells"))
  colnames(scatter_data)[c(4,5)] = c("Feat_A", "Feat_B")

  p4_plots = feat_long %>%
    split(.$Feature) %>%
    lapply(function(df) {
      ggplot(df, aes_string(
        x = colnames(feat_long)[2],
        y = colnames(feat_long)[1],
        color = "Log_Expression")) +
        geom_point(size = 2) +
        scale_color_viridis_c() +
        theme_pubr(legend = "left") +
        theme(
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)
        ) +
        ggtitle(
          paste0(unique(df$Feature)),
          subtitle = paste0("Spearman coef = ", round(coef, 2),
                            ", N cells = ", length(reg_corr_cells))
        )
    })
  p4 = ggarrange(plotlist = p4_plots)

  p5 = ggpubr::ggscatter(scatter_data, x = "Feat_A", y = "Feat_B") +
    xlab(sub(".*rna_","",rna_feat)) +
    ylab(sub(".*metabo_","",metabo_feat)) +
    ggtitle(paste0("Spearman coef = ", round(coef,2)),
            subtitle = paste0("N cells = ", length(reg_corr_cells)))  +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

  if(!is.null(cell_type)){
    final_plot_list = list(p1,p2,p3,p4,p5)
  }
  else{
    final_plot_list = list(p1,p2,p4,p5)
  }

  return(final_plot_list)
}


# Spatial Glue ------------------------------------------------------------

plot_modality_imp_per_grp = function(cell_metadata, cell_grp){
  plot_data = cell_metadata[,c("X",cell_grp,"alpha_RNA","alpha_metabo")] %>%
    dplyr::rename("Cell_grp" = cell_grp) %>%
    dplyr::filter(Cell_grp != "") %>%
    gather(key = "Modality", value = "Importance", -X, -Cell_grp)

  p = ggplot(plot_data, aes(x = Modality, y = Importance, fill = Modality)) +
    geom_violin(scale = "width", alpha = 0.7, trim = TRUE) +  # Split violins
    geom_hline(yintercept = c(0.25, 0.5, 0.75),
               linetype = "dotted",
               colour = "darkred") +
    facet_wrap(~ Cell_grp, scales = "free_y") +  # Facet by Group
    scale_fill_viridis_d() +  # Nice color scheme
    theme_pubr() +
    labs(x = "Modality", y = "Importance") +
    ggtitle("RNA vs Metabo Importance", subtitle = cell_grp) +
    theme(strip.text = element_text(size = 10, face = "bold"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  return(p)
}


# Compare corr to MISTY ---------------------------------------------------
plot_density_hm_misty = function(seurat_obj, misty_results, ROI, metabolite){
  seurat_subset = seurat_obj
  seurat_subset@graphs = list()

  seurat_subset = subset(seurat_subset, subset = region == ROI)

  misty_genes = misty_results[[ROI]] %>%
    dplyr::filter(ion == metabolite) %>%
    dplyr::pull(Predictor)


  expr_data = FetchData(seurat_subset,
                        vars = c(metabolite,misty_genes))

  binary_expr_data = binarize_dataframe(expr_data)

  codetection_mat = list()
  for (i in misty_genes){
    diff = binary_expr_data[[i]] - binary_expr_data[[paste0("metabo_",
                                                            gsub("[+-]", ".", metabolite))]]

    both_zero = binary_expr_data[[i]] == 0 & binary_expr_data[[paste0("metabo_",
                                                                      gsub("[+-]", ".", metabolite))]] == 0

    a = data.frame(var = diff)
    colnames(a)[1] = i
    codetection_mat[[i]] = a
  }
  codetection_mat = dplyr::bind_cols(codetection_mat) %>%
    as.matrix()

  # p = ComplexHeatmap::densityHeatmap(codetection_mat)

  return(codetection_mat)
}



# Compare DE to MISTY -----------------------------------------------------
plot_box_violin_feats <- function(seurat_obj, misty_de_comp, subset_cells) {
  # Check if features are in the Seurat object
  DefaultAssay(seurat_obj) = "Metabo"

  seurat_obj@graphs = list()

  features = misty_de_comp$Metabolite %>% unique()


  valid_features <- features[features %in% rownames(seurat_obj)]
  if (length(valid_features) == 0) {
    stop("None of the provided features are present in the Seurat object.")
  }

  # Subset cells
  if (!all(subset_cells %in% colnames(seurat_obj))) {
    stop("Some subset_cells are not present in the Seurat object.")
  }

  seurat_sub <- subset(seurat_obj, cells = subset_cells)
  expr_mat <- GetAssayData(seurat_sub, slot = "data")  # log-normalized data

  # Get data in long format
  expr_data <- FetchData(seurat_sub, vars = valid_features) %>%
    mutate(cell = rownames(.)) %>%
    pivot_longer(cols = all_of(valid_features), names_to = "feature", values_to = "expression") %>%
    dplyr::left_join(misty_de_comp, by = c("feature" = "Metabolite"))

  color_map <- c("#F9918A","#33C860","#81B0FF")
  names(color_map) = c("Both", "DE only", "MISTY Only")

  # Plot
  plot <- ggboxplot(expr_data, x = "feature", y = "expression",
                   fill = "Source", width = 0.7) +
    scale_fill_manual(values = color_map) +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("") + ylab("Log(Abundance)")


  gene_means <- Matrix::rowMeans(expr_mat)
  gene_vars <- apply(expr_mat, 1, var)
  mv_df <- data.frame(
    gene = rownames(expr_mat),
    mean = gene_means,
    variance = gene_vars,
    is_feature = rownames(expr_mat) %in% valid_features
  ) %>%
    dplyr::left_join(misty_de_comp, by = c("gene" = "Metabolite"))

  mv_plot <- ggplot(mv_df, aes(x = mean, y = variance)) +
    geom_point(aes(color = Source), alpha = 0.8, size = 3) +
    # scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
    theme_minimal() +
    labs(title = "Mean vs Variance of Gene Expression",
         x = "Mean Expression (log-normalized)",
         y = "Variance",
         color = "Selected Feature")

  return(list(violin_plot = plot,
              mean_variance_plot = mv_plot))
}
plot_features_mean_variance <- function(seurat_obj, misty_de_comp, metadata_col, highlight_value) {
  DefaultAssay(seurat_obj) = "Metabo"

  seurat_obj@graphs = list()

  features = misty_de_comp$Metabolite %>% unique()

  # Check if all features exist
  missing_features <- features[!features %in% rownames(seurat_obj)]
  if (length(missing_features) > 0) {
    stop("The following features are not found in the Seurat object: ", paste(missing_features, collapse = ", "))
  }

  # Check metadata column
  if (!(metadata_col %in% colnames(seurat_obj@meta.data))) {
    stop("Metadata column not found in the Seurat object.")
  }

  # Extract expression data and metadata
  expr_data <- FetchData(seurat_obj, vars = c(features, metadata_col))
  expr_long <- expr_data %>%
    pivot_longer(cols = all_of(features), names_to = "feature", values_to = "expression") %>%
    rename(group = !!metadata_col)

  # Compute mean and variance per group-feature combo
  summary_stats <- expr_long %>%
    group_by(feature, group) %>%
    summarise(
      mean_expr = mean(expression, na.rm = TRUE),
      var_expr = var(expression, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(highlight = ifelse(group == highlight_value,
                                     "ROI", "Other")) %>%
    dplyr::left_join(misty_de_comp, by = c("feature" = "Metabolite"))


  color_map <- c("#F9918A","#33C860","#81B0FF")
  names(color_map) = c("Both", "DE only", "MISTY Only")


  # Create a unique lookup of Group -> StripColorSource
  source_colors <- summary_stats %>%
    dplyr::distinct(feature, Source) %>%
    dplyr::mutate(color = color_map[Source]) %>%
    dplyr::select(feature, color) %>%
    deframe()  # makes it a named vector: Group = color

  summary_stats$feature = factor(summary_stats$feature,
                                 levels = names(source_colors))

  # Plot mean vs variance with facets
  p <- ggplot(summary_stats, aes(x = mean_expr, y = var_expr, label = group)) +
    geom_point(aes(color = highlight), size = 3, alpha = 0.5) +
    scale_color_manual(values = c("ROI" = "red", "Other" = "gray")) +
    # geom_text(vjust = -0.8, size = 3) +
    ggh4x::facet_wrap2(
      ~ feature,
      strip = ggh4x::strip_themed(
        background_x = lapply(source_colors,
                              function(fill){
                                element_rect(fill = fill, color = "black")
                              })
      ),scales = "free"
    ) +
    labs(
      title = paste("Mean vs Variance by", metadata_col),
      x = "Mean Expression",
      y = "Variance of Expression"
    ) +
    theme_pubr()

  return(p)
}
