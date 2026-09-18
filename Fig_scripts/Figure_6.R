library(tidyverse)
library(Seurat)
library(ggpubr)
library(uwot)
library(sf)

source("../Utils/utils.R")
source("../Utils/Plotting.R")
source("../Utils/colors.R")
source("../Utils/Figure_analysis.R")

main_plot_outdir = "../../../../Manuscript/Figures/Main/Figure 6"
supp_plot_outdir = "../../../../Manuscript/Figures/Supp/Figure 6"

# Load Data ---------------------------------------------------------------

metabo_RNA_seurat = readRDS("../Data/Figure 1/updated_Metabo_RNA_combined_seurat.rds")
all_sg_res = readRDS("../Data/Figure 6/updated_SG_eval_optim_res.rds")
atlas_geojson = readRDS("../Data/Figure 1/updated_atlas_geojson.rds")
joint_embed = read.csv("../Data/Figure 6/updated_Joint_embedding_best_param.csv")
brain_ontology = read.csv("../Data/Figure 3/Brain_ontology_df.csv")


# Add parent regions ------------------------------------------------------
metabo_RNA_seurat = add_parent_regions(metabo_RNA_seurat)



# Subset data for selected param comb -------------------------------------
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

joint_embed_mat = joint_embed
joint_embed_mat$X = sel_sg_res$SG_clusts_metadata$X
colnames(joint_embed_mat) = c("Cell", paste0("SG", c(1:50)))

joint_embed_mat = joint_embed_mat %>% column_to_rownames("Cell")


# UMAP joint embedding ----------------------------------------------------
library(tidytext)
SG_umap <- umap(joint_embed_mat, n_neighbors = 30,
                 min_dist = 0.2, spread = 2, metric = "euclidean")

umap_df <- as.data.frame(SG_umap)
colnames(umap_df) <- c("UMAP1", "UMAP2")

cell_region_parent = metabo_RNA_seurat@meta.data %>%
  rownames_to_column("Cell") %>%
  dplyr::select(Cell, region_parent) %>%
  dplyr::distinct()
umap_df = umap_df %>%
  rownames_to_column("Cell") %>%
  dplyr::left_join(sel_sg_res$SG_clusts_metadata, by = c("Cell" = "X")) %>%
  dplyr::mutate(SG_leiden = paste0("SG_", SG_leiden)) %>%
  dplyr::left_join(cell_region_parent) %>%
  dplyr::filter(Cell_type_coarse != "",region_parent != "Uncertain")

p2 = ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = Cell_type_coarse)) +
  geom_point(size = 1) +
  theme_pubr(legend = "right") +
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))
p2 = p2 %>% clean_labels() +
  scale_color_manual(values = assign_colors(labels = unique(umap_df$Cell_type_coarse),
                                            cols = my_cols$coarse_cell_types)) +
  NoAxes()
umap_df_with_region_names = umap_df %>%
  dplyr::left_join(brain_ontology[,c("acronym", "name")],
                   by = c("region_parent" = "acronym"))
p2.5 = ggplot(umap_df_with_region_names, aes(x = UMAP1, y = UMAP2, colour = name)) +
  geom_point(size = 1) +
  theme_pubr(legend = "right") +
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 18))
p2.5 = p2.5 %>% clean_labels() +
  scale_color_manual(values = my_cols$regions_parents) +
  NoAxes()

ct_region_parent_composition = compute_group_percent(
  umap_df, count_vars = c("region_parent", "Cell_type_coarse"), percent_group_var = "Cell_type_coarse"
) %>%
  dplyr::left_join(brain_ontology[,c("acronym", "name")],
                   by = c("region_parent" = "acronym"))

ct_region_parent_composition <- ct_region_parent_composition %>%
  mutate(region_parent_ordered = tidytext::reorder_within(name, Percent, Cell_type_coarse))

p3 = ggplot(ct_region_parent_composition, aes(x = Cell_type_coarse, y = Percent,
               fill = name, group = region_parent_ordered)) +
  geom_bar(stat = "identity", color = "white", position = position_stack()) +
  geom_text(
    data = ct_region_parent_composition %>% filter(Percent >= 10),
    aes(label = paste0(round(Percent, 1), "%")),
    position = position_stack(vjust = 0.5),
    size = 5, color = "black"
  ) +
  scale_fill_manual(values = my_cols$regions_parents) +
  coord_flip() +
  labs(x = "", y = "Percent", fill = "Parent region") +
  theme_pubr(legend = "right") +
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        axis.text = element_text(size = 14))

umaps_CT_region_parents = p2

umaps_CT_region_parents <- add_small_umap_axes(
  p = umaps_CT_region_parents,
  reduction = "umap"
)

ggsave(file.path(main_plot_outdir, "UMAP_CT_and_region_parent.pdf"),
       umaps_CT_region_parents,width = 6, height = 4)

str_region_umap_df = umap_df %>% dplyr::filter(region_parent == "STR")
p_str_region = ggplot(str_region_umap_df, aes(x = UMAP1, y = UMAP2, colour = Region)) +
  geom_point(size = 1) +
  theme_pubr() +
  NoAxes() +
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 18))

str_by_region_combined = ggarrange(plotlist = list(p3,p_str_region), nrow = 2, ncol = 1)

ggsave(file.path(supp_plot_outdir, "UMAP_STR_by_region.pdf"),
       str_by_region_combined,width = 11, height = 7)

clust_cols <- my_cols$spatialglue_clusters

umap_df$SG_leiden = factor(umap_df$SG_leiden,
                                 levels = paste0("SG_", c(0:length(unique(umap_df$SG_leiden)))))
p1 = ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = SG_leiden)) +
  geom_point(size = 1) +
  scale_color_manual(values = clust_cols) +
  theme_pubr(legend = "right") +
  NoAxes() +
  labs(color = "SG Leiden")+
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        axis.text = element_text(size = 14))

p1 <- add_small_umap_axes(
  p = p1,
  reduction = "umap"
)

ggsave(file.path(main_plot_outdir, "UMAP_SG_leiden.pdf"),
       p1,width = 6, height = 4)


p4 = ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = alpha_metabo)) +
  geom_point(size = 2) +
  scale_color_gradient2(low = scale_cols$div$low,
                        mid = scale_cols$div$mid,
                        high = scale_cols$div$high,midpoint = 0.5) +
  theme_pubr() +
  NoAxes()

p4 <- add_small_umap_axes(
  p = p4,
  reduction = "umap"
)

p5 = ggplot(umap_df, aes(x = UMAP1, y = UMAP2, colour = alpha_RNA)) +
  geom_point(size = 2) +
  scale_color_gradient2(low = scale_cols$div$low,
                        mid = scale_cols$div$mid,
                        high = scale_cols$div$high,midpoint = 0.5) +
  theme_pubr() +
  NoAxes()

umpa_importance =ggarrange(plotlist = list(p4,p5))
export_plot(plot = umpa_importance,
            path = file.path(supp_plot_outdir, "UMAP_Importance_modality.pdf"),
            sideways = T)

# Spatial plot importance -----------------------------------------------------------
modality_imp_data = sel_sg_res$SG_clusts_metadata[,c("X","x","y",
                                                     "alpha_RNA",
                                                     "alpha_metabo",
                                                     "alpha_metabo_spatial",
                                                     "alpha_metabo_feat",
                                                     "alpha_RNA_spatial",
                                                     "alpha_RNA_feat")] %>%
  gather(key = "Modality", value = "Importance", -X, -x, -y)

modality_imp_data = modality_imp_data %>%
  dplyr::filter(Modality %in% c("alpha_RNA",
                                "alpha_metabo"))

modality_imp_data$Modality = str_replace_all(modality_imp_data$Modality,
                                             c("alpha_metabo" = "Metabo",
                                               "alpha_RNA" = "RNA"))

p1 = ggplot(modality_imp_data, aes(x = -x, y = y, colour = Importance)) +
  geom_point(size = 0.2) +
  theme_pubr(legend = "bottom") +
  scale_color_gradient2(
    low = scale_cols$div$low,
    mid = scale_cols$div$mid,
    high = scale_cols$div$high,
    midpoint = 0.5
  ) +
  facet_wrap(~Modality) +
  theme(text = element_text(size = 16)) +
  NoAxes() +
  coord_fixed()

p2 = ggboxplot(data = modality_imp_data,
               x = "Modality", y = "Importance",
               fill = "Modality") +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.7, fatten = 1) +
  scale_fill_manual(values = my_cols$highlight_colors) +
  xlab("") +
  theme_pubr(legend = "none") +
  scale_x_discrete(labels = c("alpha_RNA" = "Transcriptomics",
                              "alpha_metabo" = "Metabolomics")) +
  theme(text = element_text(size = 16)) +
  rotate_x_text(angle = 45)
p = ggarrange(plotlist = list(p1,p2),widths = c(2.5,1))

ggsave(file.path(main_plot_outdir, "Spatial_Importance_modality.pdf"),
       p, width = 6, height = 4)

# Spatial polygons Importance ---------------------------------------------
SG_region_imp = sel_sg_res$SG_clusts_metadata %>%
  dplyr::select(X, Region, alpha_metabo, alpha_RNA) %>%
  dplyr::filter(Region != "") %>%
  dplyr::mutate(Imp_Diff = alpha_metabo - alpha_RNA) %>%
  dplyr::select(-alpha_RNA, -alpha_metabo) %>%
  gather(key = "Modality", value = "Importance_Diff", -X, -Region) %>%
  dplyr::group_by(Region) %>%
  dplyr::summarise(Region_Imp_Diff = mean(Importance_Diff)) %>%
  dplyr::mutate(Region = gsub("_", "-", Region))

SG_region_polygons = prep_atlas_polygons(atlas_geojson, metabo_RNA_seurat, brain_ontology) %>%
  dplyr::left_join(SG_region_imp, by = c("region" = "Region"))

p1 = ggplot(data = SG_region_polygons, aes(fill = Region_Imp_Diff)) +
  geom_sf(color = "black") +
  geom_sf_text(aes(label = region),
               color = "black",
               size = 4,
               check_overlap = T) +
  scale_fill_gradient2(
    low = scale_cols$div$low,
    mid = scale_cols$div$mid,
    high = scale_cols$div$high,midpoint = 0
  ) +
  theme_pubr(legend = "bottom") +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        aspect.ratio = 1,
        text = element_text(size = 16)) +
  labs(fill = "Importance difference\nMetabolite - gene")


ggsave(file.path(main_plot_outdir, "Polygon_Diff_Importance_modality.pdf"),
       p1, width = 8, height = 6)

# Plotting Importance distrib by group ------------------------------------
p1 = plot_modality_imp_per_grp(cell_metadata = sel_sg_res$SG_clusts_metadata,
                          cell_grp = "Region") +
  scale_fill_manual(values = my_cols$highlight_colors)
export_plot(plot = p1,
            path = file.path(supp_plot_outdir, "Importance_modality_each_region_violin.pdf"),
            sideways = T)

# Spatial plots SG leiden -------------------------------------------------
sg_clust_data = sel_sg_res$SG_clusts_metadata[,c("X","x","y",
                                                 "SG_leiden")] %>%
  dplyr::mutate(SG_leiden = paste0("SG_", SG_leiden))

clust_cols <- my_cols$spatialglue_clusters

sg_clust_data$SG_leiden = factor(sg_clust_data$SG_leiden,
                                 levels = paste0("SG_", c(0:length(unique(sg_clust_data$SG_leiden)))))

p1 = ggplot(sg_clust_data, aes(x = -x, y = y, colour = SG_leiden)) +
  geom_point(size = 0.25) +
  scale_color_manual(values = clust_cols) +
  theme_pubr(legend = "none") +
  theme(text = element_text(size = 16)) +
  NoAxes() +
  coord_fixed()

ggsave(file.path(main_plot_outdir, "Spatial_SG_leiden.pdf"),
       p1, width = 3, height = 3)

# Spatial SG leiden confusion region plot --------------------------------------
region_SG_prop = compute_dominant_sg_per_region(sel_sg_res$SG_clusts_metadata)

SG_region_polygons = prep_atlas_polygons(atlas_geojson, metabo_RNA_seurat, brain_ontology) %>%
  dplyr::left_join(region_SG_prop, by = c("region" = "Region")) %>%
  dplyr::filter(!is.na(prop))

p1 <- ggplot(SG_region_polygons, aes(fill = SG_leiden)) +
  geom_sf(color = "white") +
  scale_fill_manual(values = my_cols$spatialglue_clusters) +
  theme_pubr(legend = "right") +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    aspect.ratio = 1,
    text = element_text(size = 16)
  ) +
  labs(fill = "Dominant SG Cluster")
p2 <- ggplot(SG_region_polygons) +
  geom_sf(aes(fill = prop), color = "white", size = 0.1) +
  scale_fill_gradient(low = scale_cols$gradient$low,
                      high = scale_cols$gradient$high) +
  theme_pubr() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position = "right",
    aspect.ratio = 1,
    text = element_text(size = 16)
  ) +
  labs(fill = "Dominance strength")

combined_SG_region = ggarrange(plotlist = list(p1,p2),ncol = 2, nrow = 1)

ggsave(file.path(main_plot_outdir, "Polygon_SG_leiden_proportion.pdf"),
       combined_SG_region, width = 12, height = 4)

# Plotting importance per SG cluster violin box plot -------------------------------
plot_data = sel_sg_res$SG_clusts_metadata[,c("X","SG_leiden",
                                             "alpha_RNA","alpha_metabo")] %>%
  dplyr::rename("Cell_grp" = "SG_leiden") %>%
  dplyr::filter("SG_leiden" != "") %>%
  dplyr::mutate(Cell_grp = paste0("SG_", Cell_grp))
plot_data$Cell_grp = factor(plot_data$Cell_grp,
                            levels = paste0("SG_", c(0:length(unique(plot_data$Cell_grp)))))

plot_data$Diff_importance = plot_data$alpha_metabo - plot_data$alpha_RNA

clust_cols <- my_cols$spatialglue_clusters


p1 = ggpubr::ggviolin(data = plot_data,
                 x = "Cell_grp", y = "Diff_importance",
                 add = "boxplot",
                 fill = "Cell_grp",
                 xlab = "SpatialGlue Clusters",
                 ylab = "Importance Difference\nMetabo - RNA") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = clust_cols) +
  theme_pubr(legend = "none") +
  rotate_x_text(angle = 45) +
  theme(text = element_text(size = 16)) +
  xlab("")


ggsave(file.path(main_plot_outdir, "SG_leiden_boxplot_Importance.pdf"),
       p1, width = 6, height = 4)

# AMI plot compared to Region  --------------------------------------------
ami_res = compute_ami_by_region(metabo_RNA_seurat, include_sg = TRUE)

p = plot_ami_barplot(ami_res$plot_data)

ggsave(file.path(main_plot_outdir, "AMI_spatialglue_vs_rest_region.pdf"),
       p, width = 6, height = 4)

# Manuscript-text numbers (figure_6_numbers.csv, figure_6_dominant_SG_per_region.csv) are
# generated by Fig_scripts/Generate_Figure_Numbers.R.

# Supp Table 24 (SpatialGlue joint embedding) is now generated by
# Analysis_scripts/Downstream/Generate_Supp_Tables.R from sel_sg_res + joint_embed_mat directly.

