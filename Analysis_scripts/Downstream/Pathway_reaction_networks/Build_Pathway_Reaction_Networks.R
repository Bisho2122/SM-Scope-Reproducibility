library(tidyverse)
library(S2IsoMErData)
library(KEGGgraph)

source("../../../Utils/Plotting.R")
source("../../../Utils/utils.R")

# Data ---------------------------------------------------------------------
gene_pathway_bg = readRDS("../../../Data/Figure 4/gene_pathway_rampdb.rds")
data("Metabo_pathway_sf", package = "S2IsoMErData")
metabo_molec_mapping = read.csv("../../../Data/General/updated_annots_w_1_1_mapping.csv")
mouse_to_human_uniprot = readxl::read_xlsx("../../../Data/General/mouse_symbols_to_human_uniprot.xlsx")
RAMP_pathway_info = readRDS("../../../Data/Figure 4/RAMP_pathway_info.rds")
chebi_chemical_data = read.delim("../../../Data/General/chebi_chemical_data_April_25_2025.tsv")
Resolve_gene_uniprot = read.delim("../../../Data/Figure 1/Resolve_gene_uniprots.tsv")

# Download smpdb pathway --------------------------------------------------
path = "../../../Data/General"
message("Downloading metabolite information...")
utils::download.file(
  "https://www.smpdb.ca/downloads/smpdb_metabolites.csv.zip",
  destfile = file.path(path, "smpdb_metabolites.csv.zip")
)

message("Downloading pathway information...")
utils::download.file(
  "https://www.smpdb.ca/downloads/smpdb_pathways.csv.zip",
  destfile = file.path(path, "smpdb_pathways.csv.zip")
)

utils::unzip(
  zipfile = file.path(path, "smpdb_pathways.csv.zip"),
  exdir = file.path(path, "smpdb_pathways"),
  overwrite = TRUE
)

utils::unzip(
  zipfile = file.path(path, "smpdb_metabolites.csv.zip"),
  exdir = file.path(path, "smpdb_metabolites"),
  overwrite = TRUE
)

smpdb_pathways <-readr::read_csv(file.path(path, "smpdb_pathways/smpdb_pathways.csv"),
                  show_col_types = FALSE)

# Get datasets specific metabo and gene common pathways -------------------------------------
sf_without_adduct = sub("[+-].*", "", metabo_molec_mapping$ion) %>% unique()
metabo_pathway_bg = lapply(Metabo_pathway_sf, function(x){
  x[x %in% sf_without_adduct]
})
metabo_pathway_bg = metabo_pathway_bg[lengths(metabo_pathway_bg) > 0]

common_pathways = intersect(names(metabo_pathway_bg),
                            names(gene_pathway_bg))

to_remove_smpdb = smpdb_pathways %>%
  dplyr::filter(Subject %in% c("Disease", "Drug Action",
                                "Drug Metabolism"))

common_pathways = common_pathways[!toupper(common_pathways) %in%
                                    toupper(to_remove_smpdb$Name)]

hmdb_common_pathways = intersect(common_pathways, smpdb_pathways$Name)

# Parse Biopax from smpdb to igraph objects -------------------------------
hmdb_pathways_info = smpdb_pathways[smpdb_pathways$Name %in% hmdb_common_pathways,]
# Local cache of SMPDB biopax .owl files -- falls back to downloading from smpdb.ca if missing
local_biopax_paths = file.path(path, "smpdb_biopax")
hmdb_sif_list = list()
for (i in hmdb_common_pathways){
  if (i %in% names(hmdb_sif_list)){
    next()
  }
  else{
    PW_id = smpdb_pathways$`PW ID`[smpdb_pathways$Name == i]
    if (file.exists(file.path(local_biopax_paths, paste0(PW_id, ".owl")))){
      pax_file = paxtoolsr::readBiopax(file.path(local_biopax_paths, paste0(PW_id, ".owl")))
      sif_list = paxtoolsr::toSifnx(pax_file)
    }
    else{
      pathway_id = smpdb_pathways$`SMPDB ID`[smpdb_pathways$Name == i]
      url = paste0("https://smpdb.ca/view/", pathway_id, "/download?type=owl_markup")

      tmp = tempfile()
      download_result = tryCatch(
        download.file(url = url, destfile = tmp),
        error = function(e){NA}
      )
      if (is.na(download_result)){
        next()
      }

      pax_file = paxtoolsr::readBiopax(tmp)
      sif_list = paxtoolsr::toSifnx(pax_file)
    }
  }

  hmdb_sif_list[[i]] = sif_list
}

hmdb_rxn_nets = lapply(hmdb_sif_list, function(x){
  net = x[["edges"]] %>%
    dplyr::filter(INTERACTION_TYPE == "catalysis-precedes") %>%
    dplyr::select(PARTICIPANT_A, PARTICIPANT_B)
  colnames(net) = c("from", "to")
  net
})
names(hmdb_rxn_nets) = names(hmdb_sif_list)


# Map metabolites and gene uniprots hmdb --------------------------------------
local_metabolite_paths = file.path(path, "smpdb_metabolites")
hmdb_metabo_names_sf = list()
for (i in hmdb_common_pathways){
  pathway_id = smpdb_pathways$`SMPDB ID`[smpdb_pathways$Name == i]
  metabo_mapping_pathway = read.csv(file.path(local_metabolite_paths, paste0(pathway_id, "_metabolites.csv")))

  a = metabo_mapping_pathway$Formula
  names(a) = metabo_mapping_pathway$Metabolite.Name

  hmdb_metabo_names_sf[[i]] = a
}

hmdb_gene_names = mouse_to_human_uniprot$From
names(hmdb_gene_names) = mouse_to_human_uniprot$Entry

hmdb_gene_metabo_sif = hmdb_sif_list
for (i in hmdb_common_pathways){
 hmdb_gene_metabo_sif[[i]][["nodes"]]$PARTICIPANT = str_replace_all(hmdb_gene_metabo_sif[[i]][["nodes"]]$PARTICIPANT,
                                                             hmdb_metabo_names_sf[[i]])

 hmdb_gene_metabo_sif[[i]][["edges"]]$PARTICIPANT_A = str_replace_all(hmdb_gene_metabo_sif[[i]][["edges"]]$PARTICIPANT_A,
                                                             hmdb_metabo_names_sf[[i]])
 hmdb_gene_metabo_sif[[i]][["edges"]]$PARTICIPANT_B = str_replace_all(hmdb_gene_metabo_sif[[i]][["edges"]]$PARTICIPANT_B,
                                                               hmdb_metabo_names_sf[[i]])
}

hmdb_gene_metabo_sif = lapply(hmdb_gene_metabo_sif, function(x){
  net = x[["edges"]] %>%
    dplyr::filter(INTERACTION_TYPE %in%
                    c("controls-production-of",
                      "consumption-controlled-by")) %>%
    dplyr::filter(PARTICIPANT_A %in% sf_without_adduct | PARTICIPANT_B %in% sf_without_adduct) %>%
    dplyr::select(PARTICIPANT_A,INTERACTION_TYPE,PARTICIPANT_B)
  net
})

hmdb_gene_pathway = lapply(hmdb_sif_list, function(x){
  x[["nodes"]] %>%
    dplyr::filter(PARTICIPANT_TYPE == "ProteinReference",
                  PARTICIPANT %in% mouse_to_human_uniprot$Entry) %>%
    dplyr::left_join(mouse_to_human_uniprot[,c("From", "Entry")],
                     by = c("PARTICIPANT" = "Entry")) %>%
    dplyr::rename("mmu_gene" = "From")
})

save(hmdb_sif_list,hmdb_pathways_info,
     hmdb_rxn_nets,hmdb_gene_metabo_sif,
     hmdb_gene_pathway,
     hmdb_common_pathways, file = "../../../Data/HMDB_networks.RData")

# Download Reactome reactions chebi and uniprot -------------------------------------------------------

tmp = tempfile()
utils::download.file(
  "https://reactome.org/download/current/ChEBI2ReactomeReactions.txt",
  destfile = tmp
)
chebi2reactome_rxns = read.delim(tmp, header = F)
colnames(chebi2reactome_rxns) = c("chebi_id", "rxn_id",
                             "rxn_url", "rxn_desc", "rxn_source",
                             "species")

tmp = tempfile()
utils::download.file(
  "https://reactome.org/download/current/UniProt2ReactomeReactions.txt",
  destfile = tmp
)
uniprot2reactome_rxns = read.delim(tmp, header = F)
colnames(uniprot2reactome_rxns) = c("uniprot", "rxn_id",
                                  "rxn_url", "rxn_desc", "rxn_source",
                                  "species")


# Get Reactome dataset specific -------------------------------------------
ds_chebi_data = chebi_chemical_data[chebi_chemical_data$formula %in% sf_without_adduct,]
ds_chebi_reactome_rxns = ds_chebi_data %>%
  dplyr::filter(compound_id %in% chebi2reactome_rxns$chebi_id) %>%
  dplyr::left_join(chebi2reactome_rxns, by = c("compound_id" = "chebi_id")) %>%
  dplyr::filter(species == "Mus musculus") %>%
  dplyr::select(compound_id, formula, rxn_id, rxn_desc) %>%
  dplyr::distinct()

ds_uniprot_reactom = uniprot2reactome_rxns %>%
  dplyr::filter(uniprot %in% Resolve_gene_uniprot$Entry) %>%
  dplyr::filter(species == "Mus musculus") %>%
  dplyr::select(uniprot, rxn_id, rxn_desc)

ds_rxns = unique(ds_chebi_reactome_rxns$rxn_id)
reactome_rxns_to_pathway = list()
pb = txtProgressBar(min = 0, max = length(ds_rxns), style = 3)
counter = 0
for (i in ds_rxns){
  counter = counter + 1
  reactome_rxns_to_pathway[[i]] = Get_pathways_of_rxns(rxn_id = i)
  setTxtProgressBar(pb, counter)
}

reactome_rxns_to_pathway = reactome_rxns_to_pathway %>% dplyr::bind_rows()

pathways_of_interest = intersect(reactome_rxns_to_pathway$pathway_name,
                                 common_pathways)

reactome_pathway_to_rxns = list()
pb = txtProgressBar(min = 0, max = length(pathways_of_interest), style = 3)
counter = 0
for (i in pathways_of_interest){
  counter = counter + 1
  path_id = reactome_rxns_to_pathway$pathway_id[reactome_rxns_to_pathway$pathway_name == i] %>% unique()
  reactome_pathway_to_rxns[[i]] = Get_rxns_in_pathway(pathway_id = path_id)
  setTxtProgressBar(pb, counter)
}

reactome_pathway_to_rxns = reactome_pathway_to_rxns %>%
  dplyr::bind_rows(.id = "pathway_name")

rxn_rxn_networks = list()
pb = txtProgressBar(min = 0, max = length(pathways_of_interest), style = 3)
counter = 0
for (i in pathways_of_interest){
  counter = counter + 1
  path_id = reactome_rxns_to_pathway$pathway_id[reactome_rxns_to_pathway$pathway_name == i] %>% unique()

  rns_in_pathway = reactome_pathway_to_rxns$rxn_id[reactome_pathway_to_rxns$pathway_id == path_id] %>% unique()
  rxns_net_pathway = list()
  for (r in rns_in_pathway){
    rxns_net_pathway[[r]] = Get_rxns_preceding(rxn_id = r)
  }
  rxn_rxn_networks[[i]] = rxns_net_pathway %>% dplyr::bind_rows()

  setTxtProgressBar(pb, counter)
}



saveRDS(rxn_rxn_networks, "../../../Data/Figure 4/Reactome_rxn_rxn_networks.rds")
save(rxn_rxn_networks, reactome_pathway_to_rxns,
     pathways_of_interest,reactome_rxns_to_pathway,
     ds_chebi_data,ds_chebi_reactome_rxns,ds_uniprot_reactom, file = "../../../Data/Reactome_networks.RData")

# Download KEGG data for reactions ----------------------------------------
tmp = tempfile()
utils::download.file(
  "https://rest.kegg.jp/list/reaction",
  destfile = tmp)
kegg_rxns_desc = read.delim(tmp, header = F)
colnames(kegg_rxns_desc) = c("kegg_rxn_id", "kegg_rxn_desc")

tmp = tempfile()
utils::download.file(
  "https://rest.kegg.jp/conv/compound/chebi",
  destfile = tmp)
kegg_compound_chebi = read.delim(tmp, header = F)
colnames(kegg_compound_chebi) = c("chebi_id","kegg_cpd")

tmp = tempfile()
utils::download.file(
  "https://rest.kegg.jp/link/enzyme/mmu",
  destfile = tmp)
kegg_mmu_enzyme = read.delim(tmp, header = F)
colnames(kegg_mmu_enzyme) = c("kegg_mmu_id", "kegg_ec")

tmp = tempfile()
utils::download.file(
  "https://rest.kegg.jp/conv/mmu/uniprot",
  destfile = tmp)
kegg_mmu_uniprot = read.delim(tmp, header = F)
colnames(kegg_mmu_uniprot) = c("kegg_uniprot","kegg_mmu_id")


tmp = tempfile()
utils::download.file(
  "https://rest.kegg.jp/link/rn/ec",
  destfile = tmp)
kegg_enzyme_rxns = read.delim(tmp, header = F)
colnames(kegg_enzyme_rxns) = c("kegg_ec", "kegg_rxn_id")

tmp = tempfile()
utils::download.file(
  "https://rest.kegg.jp/link/rn/cpd",
  destfile = tmp)
kegg_compound_rxns = read.delim(tmp, header = F)
colnames(kegg_compound_rxns) = c("kegg_cpd", "kegg_rxn_id")

tmp = tempfile()
utils::download.file(
  "https://rest.kegg.jp/list/pathway/mmu",
  destfile = tmp)
kegg_mouse_pathways = read.delim(tmp, header = F)
colnames(kegg_mouse_pathways) = c("kegg_pathway_id", "kegg_pathway_name")


# Build reaction networks KEGG dataset specific --------------------------------------------
kegg_POIs = RAMP_pathway_info %>%
  dplyr::filter(pathwayName %in% common_pathways,
                type == "kegg") %>%
  dplyr::mutate(sourceId = gsub("map", "rn", sourceId)) %>%
  dplyr::pull(sourceId)

kegg_rxn_nets = list()
for (i in kegg_POIs){
  if(i %in% names(kegg_rxn_nets)){
    next()
  }
  kegg_rxn_nets[[i]] = kegg_pathway2rxn_network(pathway_id = i)
}
saveRDS(kegg_rxn_nets, "../../../Data/KEGG_rxn_networks.rds")

ds_chebi_reactome_rxns = ds_chebi_data %>%
  dplyr::filter(compound_id %in% chebi2reactome_rxns$chebi_id) %>%
  dplyr::left_join(chebi2reactome_rxns, by = c("compound_id" = "chebi_id")) %>%
  dplyr::filter(species == "Mus musculus") %>%
  dplyr::select(compound_id, formula, rxn_id, rxn_desc) %>%
  dplyr::distinct()

kegg_ds_chebi_cpds_rxns = ds_chebi_data %>%
  dplyr::mutate(chebi_id = paste0("chebi:", compound_id)) %>%
  dplyr::left_join(kegg_compound_chebi) %>%
  dplyr::filter(!is.na(kegg_cpd)) %>%
  dplyr::select(formula, chebi_id, kegg_cpd) %>%
  dplyr::left_join(kegg_compound_rxns) %>%
  dplyr::filter(!is.na(kegg_rxn_id))

kegg_ds_genes_rxns = kegg_mmu_uniprot %>%
  dplyr::mutate(kegg_uniprot = gsub("up:","",kegg_uniprot)) %>%
  dplyr::filter(kegg_uniprot %in% Resolve_gene_uniprot$Entry) %>%
  dplyr::left_join(kegg_mmu_enzyme) %>%
  dplyr::left_join(kegg_enzyme_rxns) %>%
  dplyr::filter(!is.na(kegg_ec)) %>%
  dplyr::left_join(Resolve_gene_uniprot[,c("From", "Entry")],
                   by = c("kegg_uniprot" = "Entry")) %>%
  dplyr::rename("gene" = "From")

save(kegg_rxn_nets, kegg_ds_chebi_cpds_rxns, kegg_ds_genes_rxns,
     file = "../../../Data/KEGG_Reactions_Data.RData")
