library(dplyr)
library(tidyr)

# Set to your own cluster storage path (matches baysor_dev_seg_run's -o output dir)
cluster_base_dir = Sys.getenv("SPATIAL_OMICS_CLUSTER_DIR", "/path/to/your/cluster/storage/Spatial_omics")

seg_csv = read.csv(file.path(cluster_base_dir, "Results/Baysor_segmentation/Baysor_run_w_segmask_out_dev/segmentation.csv"))


filtered = seg_csv %>%
  dplyr::filter(is_noise == "false",
                prior_segmentation != 0,
                assignment_confidence >= 0.75,
                confidence >= 0.75)

thresholds = list("is_noise" = "false",
                  "prior_segmentation" = 'not_zero',
                  "assignment_confidence" = ">=0.75",
                  "confidence" = ">=0.75")

filtered_counts = filtered %>%
  dplyr::select(gene, prior_segmentation) %>%
  dplyr::group_by(gene, prior_segmentation) %>%
  dplyr::summarise(counts = n()) %>%
  as.data.frame()
Baysor_count_mat = tidyr::pivot_wider(data = filtered_counts,
                               names_from = prior_segmentation,
                               values_from = counts)

saveRDS(list("count_mat" = Baysor_count_mat,
             "thresholds" = thresholds),
        "../../../Data/Figure 1/Baysor_filtered_count_mat.rds")
