library(RaMP)
library(RMariaDB)
library(tidyverse)
pkg.globals <- setConnectionToRaMP(dbname="ramp",username="root",conpass="ScienceBisho@2122",host = "localhost")
con = RaMP::connectToRaMP()

"%nin%" = Negate("%in%")


# Data --------------------------------------------------------------------
analyte_info = dbReadTable(con, "analyte")
source = dbReadTable(con, "source")
pathway = dbReadTable(con, "pathway")
gene_has_pathway = dbReadTable(con, "analytehaspathway")

mouse_to_human_uniprot = readxl::read_xlsx("../../../Data/General/mouse_symbols_to_human_uniprot.xlsx")

saveRDS(pathway,"../../../Data/Figure 4/RAMP_pathway_info.rds")
# Human to mouse --------------------------------------------------
gene_ramp_ids = analyte_info %>%
  dplyr::filter(type == "gene") %>%
  dplyr::pull(rampId) %>%
  unique()

gene_has_pathway = gene_has_pathway %>%
  dplyr::filter(rampId %in% gene_ramp_ids) %>%
  dplyr::left_join(pathway)

gene_source = source %>%
  dplyr::filter(geneOrCompound == "gene",
                rampId %in% gene_has_pathway$rampId,
                IDtype == "gene_symbol") %>%
  dplyr::select(-dataSource) %>%
  distinct()

gene_source_dataset = gene_source %>%
  dplyr::filter(commonName %in% mouse_to_human_uniprot$`Gene Names (primary)`) %>%
  dplyr::select(rampId, commonName) %>%
  distinct() %>%
  dplyr::left_join(mouse_to_human_uniprot[,c("From", "Gene Names (primary)")],
                   by = c("commonName" = "Gene Names (primary)"))

gene_source_dataset = gene_source_dataset[gene_source_dataset$commonName != "NANOS1",]


# Gene pathway list dataset specific --------------------------------------

Gene_pathway_bg = gene_source_dataset %>%
  dplyr::left_join(gene_has_pathway) %>%
  dplyr::group_by(pathwayName) %>%
  dplyr::summarise(
    regulons = set_names(list(From), pathwayName[1]),
    .groups = "drop"
  ) %>%
  pull(regulons)

saveRDS(Gene_pathway_bg, "../../../Data/Figure 4/gene_pathway_rampdb.rds")
