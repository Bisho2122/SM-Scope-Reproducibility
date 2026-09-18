library(tidyverse)
library(ComplexHeatmap)

# Set to your own cluster storage path (matches Retrieve_METASPACE_annotations.py's/METASPACE_Annot_QC.py's writes)
cluster_base_dir = Sys.getenv("SPATIAL_OMICS_CLUSTER_DIR", "/path/to/your/cluster/storage/Spatial_omics")

# Load data ---------------------------------------------------------------
ds_annot_res = read.csv(file.path(cluster_base_dir, "Results/Upstream/METASPACE_annotations/2022-01-11_23h33m36s.csv"))

ds_annot_res = ds_annot_res %>%
  dplyr::filter(database != "HMDB-endogenous")


# Compile QC, mainly off sample  ------------------------------------------
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

# Get annots per db -------------------------------------------------------
annots_per_db = ds_annot_res %>%
  dplyr::select(ion, database) %>%
  dplyr::distinct() %>%
  dplyr::filter(ion %in% on_sample_ions)


# Make upset plot ---------------------------------------------------------
annots_db_list = list()
for (db in unique(annots_per_db$database)){
  annots_db_list[[db]] = unique(annots_per_db$ion[annots_per_db$database == db])
}
m1 = make_comb_mat(annots_db_list, mode = "distinct")
m1 = m1[comb_size(m1) >= 5]
cs = comb_size(m1)

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

other_db_annots = annots_db_list[-which(names(annots_db_list) == "RaMP_db_brain")] %>%
  unlist() %>% unique()

db_of_interest =  annots_db_list[which(names(annots_db_list) == "RaMP_db_brain")] %>%
  unlist() %>% unique()

db_of_interest_unique = setdiff(db_of_interest, other_db_annots)

db_of_interest_unique
