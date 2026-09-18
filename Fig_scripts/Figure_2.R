library(Seurat)
library(ggpubr)
library(tidyverse)
library(pals)

source("../Utils/utils.R")
source("../Utils/Plotting.R")
source("../Utils/colors.R")
source("../Utils/Figure_analysis.R")

main_plot_outdir = "../../../../Manuscript/Figures/Main/Figure 2"
supp_plot_outdir = "../../../../Manuscript/Figures/Supp/Figure 2"

# Load Data ---------------------------------------------------------------
metabo_RNA_seurat = readRDS("../Data/Figure 1/updated_Metabo_RNA_combined_seurat.rds")
metabo_RNA_seurat@meta.data$celltype_coarse = sub("_.*", "", metabo_RNA_seurat@meta.data$Cell_type_fine)

DE_res = readRDS("../Data/Figure 2/updated_sc_DE_markers.rds")
ORA_res = readRDS("../Data/Figure 2/updated_ORA_metabotypes.rds")


# Replace formulas with mol names -----------------------------------------
metabolite_mapping_res = map_metabolite_names(metabo_RNA_seurat)
metabo_RNA_seurat = metabolite_mapping_res$seurat_obj
metabo_molec_mapping = metabolite_mapping_res$mapping

DE_res$Metabo[c("Celltypes", "Metabotype")] <- lapply(DE_res$Metabo[c("Celltypes", "Metabotype")],
                                                      map_sf_to_names)

# Colors and aesthetics config --------------------------------------------
n_metabo_types = length(unique(metabo_RNA_seurat$Metabotype))
n_celltypes = length(unique(metabo_RNA_seurat$Cell_type_fine))
n_celltypes_coarse = length(unique(metabo_RNA_seurat$celltype_coarse))
n_rna_leiden = length(unique(metabo_RNA_seurat$seurat_clusters))

metabotype_colors = assign_colors(labels = metabo_RNA_seurat$Metabotype %>%
                                    unique(),
                                  cols = my_cols$metabotypes)

celltype_colors = assign_colors(labels = metabo_RNA_seurat$Cell_type_fine %>%
                                    unique(),
                                  cols = my_cols$fine_cell_types)

celltype_coarse_colors = assign_colors(labels = metabo_RNA_seurat$celltype_coarse %>%
                                  unique(),
                                cols = my_cols$coarse_cell_types)

RNA_leiden_colors = assign_colors(labels = metabo_RNA_seurat$seurat_clusters %>%
                                         unique(),
                                       cols = my_cols$transc_leiden_clusters)

# UMAPs -------------------------------------------------------------------
metabo_RNA_seurat$Metabotype = factor(metabo_RNA_seurat$Metabotype,
                                      levels = paste0("Met.", c(1:n_metabo_types)))

metabo_metabotype_umap = DimPlot(metabo_RNA_seurat,
        reduction = "Metabo_umap",
        group.by = "Metabotype",
        pt.size = 1, cols = metabotype_colors) +
  ggtitle("Metabotype") +
  theme(text = element_text(size = 16),
        plot.title = element_blank()) +
  labs(color = "Metabotype") +
  NoAxes()

metabo_metabotype_umap <- add_small_umap_axes(
  p = metabo_metabotype_umap,
  reduction = "umap"
)

export_plot(plot = metabo_metabotype_umap,
            path = file.path(main_plot_outdir, "Metabo_metabotype_UMAP.pdf"),
            sideways = T)

RNA_leiden_umap = DimPlot(metabo_RNA_seurat,
                            reduction = "RNA_umap",
                            group.by = "seurat_clusters",
                            pt.size = 1, cols = RNA_leiden_colors) +
  labs(color = "RNA\nLeiden") +
  NoAxes() +
  theme(text = element_text(size = 16),
        plot.title = element_blank())

RNA_leiden_umap <- add_small_umap_axes(
  p = RNA_leiden_umap,
  reduction = "umap"
)

export_plot(plot = RNA_leiden_umap,
            path = file.path(supp_plot_outdir, "RNA_leiden_UMAP.pdf"),
            sideways = T)

RNA_celltype_umap = DimPlot(metabo_RNA_seurat,
        reduction = "RNA_umap",
        group.by = "Cell_type_fine",
        pt.size = 1, cols = celltype_colors) +
  theme(text = element_text(size = 16),
        plot.title = element_blank()) +
  labs(color = "Cell type") +
  NoAxes()

RNA_celltype_umap <- add_small_umap_axes(
  p = RNA_celltype_umap,
  reduction = "umap"
)

export_plot(plot = RNA_celltype_umap,
            path = file.path(main_plot_outdir, "RNA_celltype_fine_UMAP.pdf"),
            sideways = T)

RNA_celltype_coarse_umap = DimPlot(metabo_RNA_seurat,
                            reduction = "RNA_umap",
                            group.by = "celltype_coarse",
                            pt.size = 1, cols = celltype_coarse_colors) +
  NoAxes() +
  theme(text = element_text(size = 16),
        plot.title = element_blank())

RNA_celltype_coarse_umap <- add_small_umap_axes(
  p = RNA_celltype_coarse_umap,
  reduction = "umap"
)

export_plot(plot = RNA_celltype_coarse_umap,
            path = file.path(supp_plot_outdir, "RNA_celltype_coarse_UMAP.pdf"),
            sideways = T)

RNA_metabotype_umap = DimPlot(metabo_RNA_seurat,
                            reduction = "RNA_umap",
                            group.by = "Metabotype",
                            pt.size = 1, cols = metabotype_colors) +
  NoAxes() +
  labs(color = "Metabotype") +
  theme(text = element_text(size = 16),
        plot.title = element_blank())

RNA_metabotype_umap <- add_small_umap_axes(
  p = RNA_metabotype_umap,
  reduction = "umap"
)

export_plot(plot = RNA_metabotype_umap,
            path = file.path(supp_plot_outdir, "RNA_metabotype_UMAP.pdf"),
            sideways = T)

metabo_celltype_umap = DimPlot(metabo_RNA_seurat,
                              reduction = "Metabo_umap",
                              group.by = "Cell_type_fine",
                              pt.size = 1, cols = celltype_colors) +
  NoAxes() +
  theme(text = element_text(size = 16),
        plot.title = element_blank()) +
  labs(color = "Cell type")

metabo_celltype_umap <- add_small_umap_axes(
  p = metabo_celltype_umap,
  reduction = "umap"
)
export_plot(plot = metabo_celltype_umap,
            path = file.path(supp_plot_outdir, "Metabo_celltype_UMAP.pdf"),
            sideways = T)


# Spatial plots -----------------------------------------------------------
spatial_celltype = SpatialPlot(metabo_RNA_seurat,
            cols = celltype_colors,
            alpha = 1,stroke = NA,
            group.by = "Cell_type_fine",pt.size.factor = 2) +
  theme(text = element_text(size = 16),
        plot.title = element_blank()) +
  NoLegend() +
  coord_fixed()

export_plot(plot = spatial_celltype,
            path = file.path(main_plot_outdir, "RNA_celltype_spatial.pdf"),
            sideways = T)

spatial_metabotype = SpatialPlot(metabo_RNA_seurat,
            cols = metabotype_colors,
            alpha = 1,stroke = NA,
            group.by = "Metabotype",image.alpha = 0,pt.size.factor = 2) +
  theme(text = element_text(size = 16),
        plot.title = element_blank()) +
  NoLegend() +
  coord_fixed()
export_plot(plot = spatial_metabotype,
            path = file.path(main_plot_outdir, "Metabo_metabotype_spatial.pdf"),
            sideways = T)

# Stacked bar plot composition --------------------------------------------
leiden_celltype_composition = compute_group_percent(
  metabo_RNA_seurat@meta.data %>%
    dplyr::mutate(cell_type = sub("_.*", "", Cell_type_fine)),
  count_vars = c("seurat_clusters", "cell_type"),
  percent_group_var = "cell_type"
)

p1 = ggbarplot(leiden_celltype_composition,
          x = "cell_type",
          y = "Percent",
          fill = "seurat_clusters",        # Stack by 'Category'
          color = "white",          # Optional: outline bars
          palette = RNA_leiden_colors,
          position = position_stack(),  # Stack the bars
          x.text.angle = 0,          # Optional: angle x-axis text
) +
  coord_flip() +
  xlab("") +
  labs(fill = "RNA\nLeiden")

export_plot(plot = p1,
            path = file.path(supp_plot_outdir,"stack_celltype_leiden.pdf"),
            sideways = F)

# Confusion heatmap/sankey metabotypes and cell types ----------------------------
confusion_metadata = metabo_RNA_seurat@meta.data %>%
  rownames_to_column("Cell") %>%
  dplyr::select(Cell, Cell_type_fine, Metabotype)

conf_mat <- table(confusion_metadata$Cell_type_fine,
                  confusion_metadata$Metabotype)

conf_mat_norm <- prop.table(conf_mat, margin = 2) * 100

conf_df <- as.data.frame(conf_mat_norm) %>%
  dplyr::rename(Cell_type_fine = Var1, Metabotype = Var2, Percentage = Freq)

chord_mat = conf_mat %>% as.data.frame()
colnames(chord_mat) = c("Cell_type_fine", "Metabotype", "n_cells")

chord_mat = chord_mat %>%
  dplyr::left_join(conf_df) %>%
  dplyr::filter(Percentage > 0) %>%
  dplyr::select(-Percentage) %>%
  dplyr::mutate(Cell_type_fine = gsub(" ", "\n", Cell_type_fine),
                Cell_type_fine = gsub("_", "\n", Cell_type_fine)) %>%
  pivot_wider(names_from = Metabotype, values_from = n_cells, values_fill = 0) %>%
  column_to_rownames("Cell_type_fine") %>%
  as.matrix() %>%
  t()

export_graph_start(file = file.path(main_plot_outdir,"Chord_celltype_metabotype.pdf"),
                   sideways = T)
chord_plot = circlize::chordDiagram(chord_mat,transparency = 0.3,
                       directional = 1,
                       link.target.prop = T, link.sort = T,
                       scale = T,grid.col = c(metabotype_colors,
                                              celltype_colors))
export_graph_end()

# DE plots ----------------------------------------------------------------

metabo_celltype_de = filter_de_markers(DE_res$Metabo$Celltypes)
metabo_celltype_markers = de_dotplot(metabo_RNA_seurat, metabo_celltype_de,
                                     assay = "Metabo", group_by = "Cell_type_fine")

ggsave(file.path(main_plot_outdir, "Metabo_celltype_markers.pdf"),
       metabo_celltype_markers, width = 7, height = 6)

rna_metabotype_de = filter_de_markers(DE_res$RNA$Metabotype)
RNA_metabotype_markers = de_dotplot(metabo_RNA_seurat, rna_metabotype_de,
                                    assay = "RNA", group_by = "Metabotype")

ggsave(file.path(main_plot_outdir, "RNA_metabotype_markers.pdf"),
       RNA_metabotype_markers, width = 5, height = 6)


metabo_metabotype_de = filter_de_markers(DE_res$Metabo$Metabotype)
metabo_metabotype_markers = de_dotplot(metabo_RNA_seurat, metabo_metabotype_de,
                                       assay = "Metabo", group_by = "Metabotype")

ggsave(file.path(supp_plot_outdir, "Metabo_metabotype_markers.pdf"),
       metabo_metabotype_markers, width = 9, height = 9)


rna_celltype_de = filter_de_markers(DE_res$RNA$Celltypes)
RNA_celltype_markers = de_dotplot(metabo_RNA_seurat, rna_celltype_de,
                                  assay = "RNA", group_by = "Cell_type_fine", rotate_x = TRUE)

ggsave(file.path(supp_plot_outdir, "RNA_celltype_markers.pdf"),
       RNA_celltype_markers, width = 5, height = 8)



# ORA plot ----------------------------------------------------------------

p1 = S2IsoMEr::dotplot_ORA(ORA_res = ORA_res) +
  scale_color_gradient(low = scale_cols$gradient$low,
                       high = scale_cols$gradient$high) +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  xlab("")

ggsave(file.path(supp_plot_outdir, "ORA_subclass_metabotype_DE_dotplot.pdf"),
       p1, width = 8, height = 4)

# Manuscript-text numbers (figure_2_numbers.csv) are generated by
# Fig_scripts/Generate_Figure_Numbers.R.
