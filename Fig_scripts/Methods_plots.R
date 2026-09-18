library(tidyverse)
library(ComplexHeatmap)
library(ggpubr)
library(tidyheatmaps)

# Supplementary methods figures assembled from Analysis_scripts outputs

source("../Utils/utils.R")
source("../Utils/Plotting.R")
source("../Utils/colors.R")

supp_plot_outdir = "../../../../Manuscript/Figures/Supp/Methods"
dir.create(supp_plot_outdir, showWarnings = FALSE, recursive = TRUE)

# Set to your own cluster storage path (matches Retrieve_METASPACE_annotations.py's/METASPACE_Annot_QC.py's writes)
cluster_base_dir = Sys.getenv("SPATIAL_OMICS_CLUSTER_DIR", "/path/to/your/cluster/storage/Spatial_omics")

# METASPACE annotation database overlap (upset plot) ------------------------
# Ported from Analysis_scripts/Upstream/METASPACE_annotations/Compare_annots_refmet_other_dbs.R

ds_annot_res = read.csv(file.path(cluster_base_dir, "Results/Upstream/METASPACE_annotations/2022-01-11_23h33m36s.csv"))

ds_annot_res = ds_annot_res %>%
  dplyr::filter(database != "HMDB-endogenous")

QC_res_path = file.path(cluster_base_dir, "Results/Upstream/METASPACE_annotations/Annotation_QC")
input_files = list.files(QC_res_path, pattern = "Annot_QC_2022")

all_off_sample = list()
for (i in input_files){
  data = read.csv(file.path(QC_res_path,i))
  db = sub(".csv.*", "",sub(".*Annot_QC_2022-01-11_23h33m36s_", "", i))
  all_off_sample[[db]] = data
}
all_off_sample = all_off_sample %>% dplyr::bind_rows()

on_sample_ions = all_off_sample %>%
  dplyr::filter(off_sample_label != "off") %>%
  dplyr::mutate(ion = paste0(formula, adduct, "-")) %>%
  dplyr::pull(ion) %>%
  unique()

annots_per_db = ds_annot_res %>%
  dplyr::select(ion, database) %>%
  dplyr::distinct() %>%
  dplyr::filter(ion %in% on_sample_ions)

annots_db_list = list()
for (db in unique(annots_per_db$database)){
  annots_db_list[[db]] = unique(annots_per_db$ion[annots_per_db$database == db])
}
m1 = make_comb_mat(annots_db_list, mode = "distinct")
m1 = m1[comb_size(m1) >= 5]
cs = comb_size(m1)

export_graph_start(file = file.path(supp_plot_outdir, "METASPACE_db_annotation_upset.pdf"),
                   sideways = T)
ht = UpSet(m1, comb_order = order(comb_size(m1), decreasing = T),
      top_annotation = HeatmapAnnotation(
        "Intersections" = anno_barplot(cs,
                                             ylim = c(0, max(cs)*1.1),
                                             border = FALSE,
                                             gp = gpar(fill = "black"),
                                             height = unit(4, "cm")
        ),
        annotation_name_side = "left",
        annotation_name_rot = 90),
      pt_size =  unit(5, "mm"))

ht = draw(ht)
od = column_order(ht)
decorate_annotation("Intersections", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"),
            default.units = "native", just = c("left", "bottom"),
            gp = gpar(fontsize = 12, col = "#404040"), rot = 45)
})
export_graph_end()

# Allen atlas affine registration RMSE (line plot) --------------------------
# Ported from Analysis_scripts/Upstream/Allen_Atlas_to_cells/plot_RMSE_affine_atlas.R

RMSE_Res = read.csv("../Data/Figure 1/Atlas_sections_RMSE_results.csv")

RMSE_Res = RMSE_Res[order(RMSE_Res$Section),]

RMSE_Res$Section = gsub("Section_", "Plane.",RMSE_Res$Section)
RMSE_Res$section_num = sub(".*Plane.","",RMSE_Res$Section) %>% as.numeric()
RMSE_Res$section_num = RMSE_Res$section_num + 1
RMSE_Res$Section = paste0("Plane.", RMSE_Res$section_num)

atlas_rmse_plot = RMSE_Res %>%
  ggplot(aes(x = Section, y = RMSE, group = 1)) +
  geom_point(size = 3) +
  geom_line() +
  theme_pubclean() +
  rotate_x_text(45) +
  xlab("") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) +
  ylim(c(200,1000))

export_plot(plot = atlas_rmse_plot,
            path = file.path(supp_plot_outdir, "Atlas_affine_registration_RMSE.pdf"))

# SpatialGlue parameter optimization (heatmap + line plot) ------------------
# Ported from Analysis_scripts/Downstream/SpatialGlue/SpatialGlue_param_optim.R
# (plotting only; the raw per-parameter-combination sweep CSVs are aggregated
# into the Data/Figure 6/ objects read below by that source script)

all_sg_res = readRDS("../Data/Figure 6/SG_eval_res_combined.rds")
metrics_data = all_sg_res$SCIB_metrics %>%
  dplyr::select(-X) %>%
  dplyr::mutate(param_comb = sub(".csv.*", "", param_comb)) %>%
  tidyr::separate(param_comb,
                  into = c("Dim", "K_feat",
                           "metric_feat", "Loss_weights"),
                  sep = "_", remove = F) %>%
  dplyr::mutate(Dim = sub("dim", "", Dim),
                K_feat = sub("k", "", K_feat),
                metric_feat = sub("metric", "", metric_feat),
                Loss_weights = sub("wf", "", Loss_weights)) %>%
  dplyr::mutate(Metric_label = paste0(Metric, ";", Label))

export_graph_start(file = file.path(supp_plot_outdir, "SpatialGlue_param_optim_heatmap.pdf"),
                   sideways = T)
tidyheatmaps::tidyheatmap(df = metrics_data,
            rows = param_comb,
            columns = Metric_label,
            values = Value,
            scale = "column",
            annotation_col = c(Label, Metric),
            annotation_row = c(Dim,K_feat, metric_feat, Loss_weights),
            gaps_col = Label,
            gaps_row = Dim,
            show_rownames = F, show_colnames = F,
            border_color = NA,display_numbers = T,
            fontsize_number = 12, cluster_rows = T,
            colors = c("#145afc","#ffffff","#ee4445"),
            cutree_rows = 18,fontsize = 12)
export_graph_end()

scaled_metrics = metrics_data %>%
  dplyr::group_by(Metric_label) %>%
  dplyr::mutate(metric_score = as.numeric(scale(Value))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(param_comb) %>%
  dplyr::mutate(n_positive_metrics = length(which(metric_score > 0))) %>%
  dplyr::ungroup()

filt_all_pos = scaled_metrics %>%
  dplyr::filter(n_positive_metrics == length(unique(scaled_metrics$Metric_label)))

filt_avg_labels = filt_all_pos %>%
  dplyr::group_by(param_comb, Metric) %>%
  dplyr::mutate(mean_metric_score = mean(Value)) %>%
  dplyr::select(-Label, -Metric_label, -Value, -metric_score) %>%
  dplyr::distinct() %>%
  dplyr::ungroup()

spatialglue_param_optim_lineplot = filt_avg_labels %>%
  ggplot(aes(x = param_comb, y = mean_metric_score, colour = Metric,
             group = Metric)) +
  geom_point(size = 3) +
  geom_line() +
  coord_flip() +
  theme_pubr() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14))

export_plot(plot = spatialglue_param_optim_lineplot,
            path = file.path(supp_plot_outdir, "SpatialGlue_param_optim_lineplot.pdf"))
