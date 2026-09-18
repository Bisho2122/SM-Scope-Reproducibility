library(tidyverse)
library(ggpubr)
require(sf)
library(Seurat)
require(DIALOGUE)
library(ggrepel)
library(uwot)

source("../Utils/utils.R")
source("../Utils/Plotting.R")
source("../Utils/colors.R")
source("../Utils/Figure_analysis.R")

main_plot_outdir = "../../../../Manuscript/Figures/Main/Figure 5"
supp_plot_outdir = "../../../../Manuscript/Figures/Supp/Figure 5"

# Load Data ---------------------------------------------------------------

metabo_RNA_seurat = readRDS("../Data/Figure 1/updated_Metabo_RNA_combined_seurat.rds")
atlas_geojson = readRDS("../Data/Figure 1/updated_atlas_geojson.rds")
brain_ontology = read.csv("../Data/Figure 3/Brain_ontology_df.csv")
dialogue_res = readRDS("../Data/Figure 5/updated_DIALOGUE_res/updated_CT_5_main_out.rds")
spatial_niches = readRDS("../Data/Figure 1/spatial_niches_regions.rds")
spatial_niches_plots = readRDS("../Data/Figure 1/spatial_niches_all_clust_results.rds")

all_DE_res = readRDS("../Data/Figure 2/updated_sc_DE_markers.rds")

metabolite_mapping_res = map_metabolite_names(metabo_RNA_seurat, rename_assay_rownames = FALSE)
metabo_molec_mapping = metabolite_mapping_res$mapping


# Add parent regions ------------------------------------------------------
metabo_RNA_seurat = add_parent_regions(metabo_RNA_seurat)


# Invert X for spatial niches ---------------------------------------------
spatial_niches_plots = lapply(spatial_niches_plots, function(x){
  x$plot$data$x = -1 * x$plot$data$x
  x
})


# Description Regions into niches plots -----------------------------------
n_niches_region = spatial_niches %>%
  dplyr::select(Region, Niche) %>%
  dplyr::distinct() %>%
  dplyr::count(Region)

p1 = n_niches_region %>%
  ggplot(aes(x = reorder(Region,n), y = n)) +
  geom_point(size = 3) +
  theme_pubr() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 16)) +
  xlab("Regions") +
  ylab("Number of Niches") +
  geom_hline(yintercept = 4, linetype = "dashed", col = "darkred",
             lwd = 0.75) +
  geom_hline(yintercept = c(3,6), linetype = "dotted", col = "blue",
             lwd = 0.75)

ggsave(file.path(supp_plot_outdir, "Number_Niches_per_region_plot.pdf"),
       p1, width = 6, height = 5)

# Spatial niches plot -----------------------------------------------------
largest_niches_plots = spatial_niches_plots[c("Isocortex", "CP", "OLF", "LHA")] %>%
  lapply(function(x){x$plot +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            text = element_text(size = 16)) +
      xlab("X") + ylab("Y")})
p1 = ggarrange(plotlist = largest_niches_plots, labels = c("Isocortex", "CP", "OLF", "LHA"),
          label.x = 0.5)

ggsave(file.path(main_plot_outdir, "Largest_spatial_niches.pdf"),
       p1, width = 8, height = 8)

# UMAP plot for MCP scores ------------------------------------------------
MCP_cell_embed = get_mcp_cell_embed(dialogue_res)

Niche_region_meta = MCP_cell_embed %>%
  dplyr::select(Niche, Region) %>%
  dplyr::distinct() %>%
  dplyr::left_join(metabo_RNA_seurat@meta.data[,c("Niche", "region_parent")]) %>%
  dplyr::filter(region_parent != "Uncertain") %>%
  dplyr::distinct()

MCP_niche_embed = MCP_cell_embed[,c("Cell", "Niche",paste0("MCP", c(1:10)))] %>%
  gather(key = "MCP", value = "MCP_score", -Cell, -Niche) %>%
  dplyr::group_by(Niche, MCP) %>%
  dplyr::summarise(MCP_score = mean(MCP_score)) %>%
  pivot_wider(values_from = MCP_score, names_from = MCP) %>%
  column_to_rownames("Niche")
MCP_niche_embed = MCP_niche_embed[,paste0("MCP", c(1:10))]

MCP_umap <- umap(MCP_niche_embed, n_neighbors = 5,
                 min_dist = 0.2, spread = 2, metric = "euclidean")

umap_df <- as.data.frame(MCP_umap)
colnames(umap_df) <- c("UMAP1", "UMAP2")

umap_df = umap_df %>%
  rownames_to_column("Niche") %>%
  dplyr::left_join(Niche_region_meta) %>%
  dplyr::left_join(brain_ontology[,c("acronym", "name")], by = c("region_parent" = "acronym"))

p1 = ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = name)) +
  geom_point(size = 1) +
  labs(colour = "Parent region")
p1 = p1 +
  scale_color_manual(values = assign_colors(labels = unique(umap_df$name),
                                            cols = my_cols$regions_parents)) +
  theme_pubr(legend = "right") +
  theme(text = element_text(size = 16),
        plot.title = element_blank()) +
  NoAxes()

p1 <- add_small_umap_axes(
  p = p1,
  reduction = "umap"
)

ggsave(file.path(main_plot_outdir, "UMAP_MCP_space_niches.pdf"),
       p1, width = 6, height = 4)

# Spatial plot MCP scores -----------------------------------------
MCP_cell_embed = get_mcp_cell_embed(dialogue_res)

MCP_region_scores = MCP_cell_embed[,c("Cell", paste0("MCP", 1:10),
                                      "Region")] %>%
  gather(key = "MCP", value = "MCP_Score", -Cell, -Region) %>%
  dplyr::group_by(Region, MCP) %>%
  dplyr::summarise(Region_MCP_Score = mean(MCP_Score)) %>%
  dplyr::mutate(Run = "New")

region_polygons = prep_atlas_polygons(atlas_geojson, metabo_RNA_seurat, brain_ontology) %>%
  dplyr::left_join(MCP_region_scores, by = c("region" = "Region"))


MCPs = unique(region_polygons$MCP)
MCPs = MCPs[!is.na(MCPs)] %>%
  factor(levels = paste0("MCP", 1:10))

MCP_spatial_plots = list()
for (i in paste0("MCP", 1:10)){
  filt = region_polygons %>%
    dplyr::filter(MCP == i)
  filt$Region_MCP_Score[abs(filt$Region_MCP_Score) < 1] = 0
  p = ggplot(data = filt, aes(fill = Region_MCP_Score)) +
    geom_sf(color = "lightgrey") +
    scale_fill_gradient2(
      low = scale_cols$div$low,
      mid = scale_cols$div$mid,
      high = scale_cols$div$high,
      midpoint = 0
    ) +
    theme_pubr(legend = "right") +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          aspect.ratio = 1,
          text = element_text(size = 16)) +
    labs(fill=i)
  MCP_spatial_plots[[i]] = p
}

combined = cowplot::plot_grid(plotlist = MCP_spatial_plots)


ggsave(file.path(main_plot_outdir, "all_MCP_scores_spatial.pdf"),
       combined, width = 12, height = 6)

# MCP - MCP correlation ---------------------------------------------------
mcp_niche_corr_mat = cor(MCP_niche_embed)
p1 = ggcorrplot::ggcorrplot(mcp_niche_corr_mat, hc.order = TRUE, type = "lower",
           outline.col = "white", lab = T,colors = scale_cols$div,
           legend.title = "Pearson \ncoefficient")

ggsave(file.path(supp_plot_outdir, "Niche_MCP_MCP_corr_plot.pdf"),
       p1, width = 6, height = 5)


# DIALOGUE package plots MCP2-------------------------------------------------
library(DIALOGUE)
dialogue_res$R$param$averaging.function = colMeans2
ct_corr_matrix = DIALOGUE:::DIALOGUE.plot.av(dialogue_res$R, MCPs = 1:10)

mcp2_corr_ct = cor(ct_corr_matrix[[2]])
colnames(mcp2_corr_ct) = c("Astrocytes", "Excitatory","Inhibitory",
                           "Microglia", "Oligodendrocytes")
rownames(mcp2_corr_ct) = c("Astrocytes", "Excitatory","Inhibitory",
                           "Microglia", "Oligodendrocytes")

p1 = ggcorrplot::ggcorrplot(mcp2_corr_ct, hc.order = TRUE, type = "lower",
                            outline.col = "white",
                            lab = T,colors = scale_cols$div,
                            legend.title = "Pearson \ncoefficient",
                            title = "MCP2",tl.cex = 16,pch.cex = 16,
                            lab_size = 5) +
  theme(text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5))

ggsave(file.path(main_plot_outdir,"MCP2_celltype_corr_pairs_plot.pdf"),
       p1, width = 6, height = 5)


# Spatial celltype plot MCP2 ----------------------------------------------
mcp2_cell_embed = MCP_cell_embed
mcp2_cell_embed$cell.type = str_replace_all(mcp2_cell_embed$cell.type, c("Excit" = "Excitatory",
                                             "Inhib" = "Inhibitory"))

p1 = ggplot(mcp2_cell_embed, aes(x = -coor.X, y = coor.Y, colour = MCP2)) +
  geom_point() +
  theme_pubr() +
  scale_color_gradient2(
    low = scale_cols$div$low,
    mid = scale_cols$div$mid,
    high = scale_cols$div$high,      # Positive values
    midpoint = 0      # Center the scale at zero
  ) +
  facet_wrap(~cell.type, scales = "free") +
  xlab("") +
  ylab("") + NoAxes()

export_plot(plot = p1,
            path = file.path(supp_plot_outdir, "MCP2_score_Celltype_spatial.pdf"),
            sideways = T)

# Spatial Region by celltype plot all MCPs ----------------------------------------------
mcp_all_cell_embed = MCP_cell_embed
mcp_all_cell_embed$cell.type = str_replace_all(mcp_all_cell_embed$cell.type, c("Excit" = "Excitatory",
                                             "Inhib" = "Inhibitory"))
cell_embed_MCP_long = mcp_all_cell_embed[,c("Cell", "cell.type", "Region","Niche",
                                    "coor.X","coor.Y",
                                    paste0("MCP", 1:10))] %>%
  gather(key = "Program", value = "Score", -Cell, -cell.type, -Region,
         -Niche, -coor.X, -coor.Y)
cell_embed_MCP_long$Program = factor(cell_embed_MCP_long$Program,
                                     levels = paste0("MCP", 1:10))
p1 = ggplot(cell_embed_MCP_long, aes(x = -coor.X, y = coor.Y, colour = Score)) +
  geom_point() +
  theme_pubr() +
  scale_color_gradient2(
    low = scale_cols$div$low,
    mid = scale_cols$div$mid,
    high = scale_cols$div$high,
    midpoint = 0
  ) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) +
  facet_grid(cell.type~Program, scales = "free") +
  NoAxes()

export_plot(plot = p1,
            path = file.path(supp_plot_outdir, "Celltype_by_MCP_spatial.pdf"),
            sideways = T)

# Metabolites in MCP2 up all cell types -----------------------------------
dialogue_sig_min5 = MCP_signature_to_df(sig_list = dialogue_res$R$sig2,
                                   min_celltype = 5) %>%
  dplyr::select(-cell_type) %>%
  dplyr::distinct()

MCP2_metabos = dialogue_sig_min5$metabolite[dialogue_sig_min5$MCP == "MCP2"]

mcp2_metabo_sign_mat = dialogue_sig_min5 %>%
  dplyr::mutate(val = ifelse(Direction == "down", -1, 1)) %>%
  dplyr::select(metabolite, MCP, val) %>%
  dplyr::filter(metabolite %in% MCP2_metabos) %>%
  dplyr::left_join(metabo_molec_mapping, by = c("metabolite" = "ion")) %>%
  dplyr::select(mol_name_one_to_one, val, MCP) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(names_from = MCP, values_from = val, values_fill = 0) %>%
  column_to_rownames("mol_name_one_to_one") %>%
  as.matrix()

pdf(file.path(main_plot_outdir,"MCP2_22_metabolites_across_MCPs.pdf"),
    width = 5,height = 5)
ComplexHeatmap::pheatmap(mcp2_metabo_sign_mat,color = c(scale_cols$div$low, scale_cols$div$mid,
                                   scale_cols$div$high),
                   cluster_rows = T, cluster_cols = T,treeheight_col = 8,
                   treeheight_row = 8, name = "MCP\nscore sign",
                   row_title_gp = gpar(fontsize = 15),
                   column_title_gp = gpar(fontsize = 15),border_color = "white")
dev.off()

# Venn diagram MCP 2 markers and metabotype 2 markers ---------------------

Met1_markers = all_DE_res$Metabo$Metabotype %>%
  dplyr::filter(avg_log2FC >= 1,
                p_val_adj < 0.05,
                cluster == "Met.1") %>%
  dplyr::pull(gene)

p1 = ggvenn::ggvenn(list("MCP 2" = MCP2_metabos,
                    "Metabotype 1" = Met1_markers),
               fill_color = my_cols$highlight_colors,
               stroke_color = NA,show_percentage = F,
               set_name_size = 4,text_size = 10,auto_scale = F)

ggsave(file.path(main_plot_outdir, "MCP2_metabotype2_venn.pdf"),
       p1, width = 4,height = 3)

# Confusion map MCPs and Metabotype ---------------------------------------
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

conf_mat_MCP_metabotype$MCP = factor(conf_mat_MCP_metabotype$MCP,
                                     levels = paste0("MCP", 1:10))
conf_mat_MCP_metabotype$Metabotype = factor(conf_mat_MCP_metabotype$Metabotype,
                                     levels = paste0("Met.", 1:length(unique(conf_mat_MCP_metabotype$Metabotype))))
p1 = conf_mat_MCP_metabotype %>%
  ggplot(aes(x = MCP, y = as.factor(Metabotype), fill = prop_common)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(prop_common,2)),
            color = "black",
            size = 4) +
  scale_fill_gradient(low = scale_cols$gradient$low,
                      high = scale_cols$gradient$high) +
  theme_minimal() +
  labs(x = "MCP", y = "Metabotype", fill = "Prop") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(supp_plot_outdir, "MCP_by_metabotype_cell_proportion.pdf"),
       p1,width = 6, height = 6)

# Supp Table 23 (MCP-metabotype overlap) is now generated by
# Analysis_scripts/Downstream/Generate_Supp_Tables.R, replicating this computation.

# MCP to Metabotype signature comparison (select 1 for main figure) -------------------------------------------------
dialogue_sig_min4 = MCP_signature_to_df(sig_list = dialogue_res$R$sig2,
                                   min_celltype = 4) %>%
  dplyr::select(-cell_type) %>%
  dplyr::distinct()
MCP_metaboytype_cells_mapping = MCP_cells_assign %>%
  dplyr::left_join(metabotype_meta)

all_metabotypes = unique(metabo_RNA_seurat@meta.data$Metabotype)

MCP_metabotype_sigs = list()
for (i in all_metabotypes){
  MCP_metabotype_sigs[[i]] = Get_metabotype_MCP_markers(dialogue_results = dialogue_res$full_RA_res,
                                                        MCP_metabotype_cells = MCP_metaboytype_cells_mapping,
                                                        metabotype_DE = all_DE_res$Metabo$Metabotype,
                                                        metabotype = i,
                                                        include_DE_overlap = F,
                                                        min_MCP_metabotype_prop = 0.8)
}

p1 = MCP_metabotype_sigs$Met.1$plots$spatial
p1[[1]] = p1[[1]] + scale_fill_manual(values = c(scale_cols$div$low,scale_cols$div$mid))
p1[[2]] = p1[[2]] + scale_fill_manual(values = c(scale_cols$div$low,scale_cols$div$mid))


export_plot(plot = p1,
            path = file.path(main_plot_outdir, "MCP2_metabotype2_spatial.pdf"),
            sideways = T)


p1 = MCP_metabotype_sigs$Met.1$plots$expr_boxplot
p1$data = p1$data %>%
  dplyr::left_join(metabo_molec_mapping,by = c("metabolite" = "ion")) %>%
  dplyr::select(-metabolite) %>%
  dplyr::rename("metabolite" = "mol_name_one_to_one")
p1$layers[[1]]$data = p1$layers[[1]]$data %>%
  dplyr::left_join(metabo_molec_mapping,by = c("metabolite" = "ion")) %>%
  dplyr::select(-metabolite) %>%
  dplyr::rename("metabolite" = "mol_name_one_to_one")

p1$data$Cell_type = str_replace_all(p1$data$Cell_type,
                                    c("Excit" = "Excitatory",
                                      "Inhib" = "Inhibitory"))

mcp2_not_metabotype2_boxplot = ggplot(p1$data, aes(x = Cell_type, y = exp, fill = Cell_type)) +
  geom_boxplot(width = 0.7, outlier.alpha = 0.25) +
  facet_wrap(~ metabolite, scales = "free_y") +  # each facet only that metabolite
  scale_fill_manual(values = my_cols$coarse_cell_types[c(1:3,6,5)]) +
  labs(x = NULL, y = "log10(abundance)") +
  theme_pubr(legend = "top") +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(6, "pt"),
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 16)
  )

ggsave(file.path(main_plot_outdir, "MCP2_not_metabotype2_boxplot_metabolites.pdf"),
       mcp2_not_metabotype2_boxplot, width = 14, height = 8)


p2 = MCP_metabotype_sigs$Met.1$plots$expr_boxplot_metabotype_only
p2$data = p2$data %>%
  dplyr::left_join(metabo_molec_mapping,by = c("metabolite" = "ion")) %>%
  dplyr::select(-metabolite) %>%
  dplyr::rename("metabolite" = "mol_name_one_to_one")
p2$layers[[1]]$data = p2$layers[[1]]$data %>%
  dplyr::left_join(metabo_molec_mapping,by = c("metabolite" = "ion")) %>%
  dplyr::select(-metabolite) %>%
  dplyr::rename("metabolite" = "mol_name_one_to_one")

p2$data$Cell_type = str_replace_all(p2$data$Cell_type,
                                    c("Excit" = "Excitatory",
                                      "Inhib" = "Inhibitory"))

mcp2_only_metabotype2_boxplot = ggplot(p2$data, aes(x = Cell_type, y = exp, fill = Cell_type)) +
  geom_boxplot(width = 0.7, outlier.alpha = 0.25) +
  facet_wrap(~ metabolite, scales = "free_y") +  # each facet only that metabolite
  scale_fill_manual(values = my_cols$coarse_cell_types[c(1:3,6,5)]) +
  labs(x = NULL, y = "log10(abundance)") +
  theme_pubr(legend = "none") +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(6, "pt"),
    strip.text = element_text(size = 14),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 16)
  )

ggsave(file.path(main_plot_outdir, "MCP2_only_metabotype2_boxplot_metabolites.pdf"),
       mcp2_only_metabotype2_boxplot, width = 5, height = 5)

# Manuscript-text numbers (figure_5_numbers.csv) are generated by
# Fig_scripts/Generate_Figure_Numbers.R.

