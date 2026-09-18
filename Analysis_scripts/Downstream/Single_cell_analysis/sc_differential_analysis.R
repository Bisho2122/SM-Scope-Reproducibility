library(Seurat)
library(tidyverse)


# Load Data ---------------------------------------------------------------

metabo_RNA_seurat = readRDS("../../../Data/Figure 1/updated_Metabo_RNA_combined_seurat.rds")


# Metabolites X cell types ------------------------------------------------
metabo_celltypes_DE = FindAllMarkers(object = metabo_RNA_seurat,
                                     assay = "Metabo",
                                     group.by = "Cell_type_fine",
                                     only.pos = F,
                                     logfc.threshold = 0.25,
                                     min.cells.group = 20,
                                     min.cells.feature = 20)

metabo_celltypes_DE$cluster = str_replace_all(metabo_celltypes_DE$cluster,
                                                c("Excit_CTX" = "Excitatory_VGLUT1",
                                                  "Inhib_CTX" = "Inhibitory_cortex",
                                                  "Excit_VGLUT2" = "Excitatory_VGLUT2",
                                                  "Excit_HPF" = "Excitatory_HPF",
                                                  "Inhib_Interneurons" = "Inhibitory_interneurons",
                                                  "Ventricular systems Glial cells" = "Ependymal cells")) %>%
  as.factor()

# Metabolites X Metabotype ------------------------------------------------
metabo_Metabotype_DE = FindAllMarkers(object = metabo_RNA_seurat,
                                     assay = "Metabo",
                                     group.by = "Metabotype",
                                     only.pos = F,
                                     logfc.threshold = 0.25,
                                     min.cells.group = 20,
                                     min.cells.feature = 20,
                                     min.pct = 0.1)


# Metabolites X Region ----------------------------------------------------
metabo_region_DE = FindAllMarkers(object = metabo_RNA_seurat,
                                  assay = "Metabo",
                                  group.by = "region",
                                  only.pos = F,
                                  logfc.threshold = 0.25,
                                  min.cells.group = 20,
                                  min.cells.feature = 20,
                                  min.pct = 0.1)


# Metabolites X CT_Region -------------------------------------------------
metabo_CT_region_DE = FindAllMarkers(object = metabo_RNA_seurat,
                                  assay = "Metabo",
                                  group.by = "CT_Region",
                                  only.pos = F,
                                  logfc.threshold = 0.25,
                                  min.cells.group = 20,
                                  min.cells.feature = 20,
                                  min.pct = 0.1)



# Genes X Metabotype ------------------------------------------------------
RNA_Metabotype_DE = FindAllMarkers(object = metabo_RNA_seurat,
                                      assay = "RNA",
                                      group.by = "Metabotype",
                                      only.pos = F,
                                      logfc.threshold = 0.25,
                                      min.cells.group = 20,
                                      min.cells.feature = 20,
                                      min.pct = 0.1)


# Genes X Celltype --------------------------------------------------------
RNA_celltypes_DE = FindAllMarkers(object = metabo_RNA_seurat,
                                   assay = "RNA",
                                   group.by = "Cell_type_fine",
                                   only.pos = F,
                                   logfc.threshold = 0.25,
                                   min.cells.group = 20,
                                   min.cells.feature = 20,
                                   min.pct = 0.1)

RNA_celltypes_DE$cluster = str_replace_all(RNA_celltypes_DE$cluster,
                                              c("Excit_CTX" = "Excitatory_VGLUT1",
                                                "Inhib_CTX" = "Inhibitory_cortex",
                                                "Excit_VGLUT2" = "Excitatory_VGLUT2",
                                                "Excit_HPF" = "Excitatory_HPF",
                                                "Inhib_Interneurons" = "Inhibitory_interneurons",
                                                "Ventricular systems Glial cells" = "Ependymal cells")) %>%
  as.factor()


# Genes X Region ----------------------------------------------------------
RNA_region_DE = FindAllMarkers(object = metabo_RNA_seurat,
                                  assay = "RNA",
                                  group.by = "region",
                                  only.pos = F,
                                  logfc.threshold = 0.25,
                                  min.cells.group = 20,
                                  min.cells.feature = 20,
                                  min.pct = 0.1)


# Genes X CT_Region -------------------------------------------------------
RNA_CT_region_DE = FindAllMarkers(object = metabo_RNA_seurat,
                               assay = "RNA",
                               group.by = "CT_Region",
                               only.pos = F,
                               logfc.threshold = 0.25,
                               min.cells.group = 20,
                               min.cells.feature = 20,
                               min.pct = 0.1)

# Save Marker results -----------------------------------------------------
if (!file.exists("../../../Data/Figure 2/updated_sc_DE_markers.rds")){
  all_DE_markers = list("RNA" = list("Celltypes" = RNA_celltypes_DE,
                                     "Metabotype" = RNA_Metabotype_DE,
                                     "Region" = RNA_region_DE,
                                     "CT_Region" = RNA_CT_region_DE),
                        "Metabo" = list("Celltypes" = metabo_celltypes_DE,
                                        "Metabotype" = metabo_Metabotype_DE,
                                        "Region" = metabo_region_DE,
                                        "CT_Region" = metabo_CT_region_DE))
} else{
  all_DE_markers = readRDS("../../../Data/Figure 2/updated_sc_DE_markers.rds")
  all_DE_markers[["RNA"]][["Region"]] = RNA_region_DE
  all_DE_markers[["Metabo"]][["Region"]] = metabo_region_DE

  all_DE_markers[["RNA"]][["CT_Region"]] = RNA_CT_region_DE
  all_DE_markers[["Metabo"]][["CT_Region"]] = metabo_CT_region_DE
}


saveRDS(all_DE_markers, "../../../Data/Figure 2/updated_sc_DE_markers.rds")

# Supp Tables 10 and 13 (metabotype/region metabolite markers) are now generated by
# Analysis_scripts/Downstream/Generate_Supp_Tables.R from this script's saved Data/ object.

