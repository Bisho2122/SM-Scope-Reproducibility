library(tidyverse)
library(igraph)

source("../../../Utils/utils.R")

#Load Data --------------------------------------------------------------------
load("../../../Data/Reactome_networks.RData")
load("../../../Data/HMDB_networks.RData")
load("../../../Data/KEGG_Reactions_Data.RData")
Resolve_gene_uniprot = read.delim("../../../Data/Figure 1/Resolve_gene_uniprots.tsv")
metabo_molec_mapping = read.csv("../../../Data/General/updated_annots_w_1_1_mapping.csv")

# Get metabo-rxns and gene-rxns for each pathway ----------------------------------
sf_without_adduct = sub("[+-].*", "", metabo_molec_mapping$ion) %>% unique()

reactome_gene_metabo = list()
for (i in 1:length(rxn_rxn_networks)){
  pathway_nodes = rxn_rxn_networks[[i]] %>% unlist() %>% as.character()
  pathway_genes_rxns = ds_uniprot_reactom %>%
    dplyr::filter(rxn_id %in% pathway_nodes) %>%
    dplyr::left_join(Resolve_gene_uniprot[,c("From", "Entry")],
                     by = c("uniprot" = "Entry")) %>%
    dplyr::filter(!is.na(From)) %>%
    dplyr::rename("gene" = "From")

  pathway_metabo_rxns = ds_chebi_reactome_rxns %>%
    dplyr::filter(rxn_id %in% pathway_nodes,
                  formula %in% sf_without_adduct) %>%
    dplyr::rename("metabo" = "formula")

  if(nrow(pathway_metabo_rxns) == 0 || nrow(pathway_genes_rxns) == 0){
    reactome_gene_metabo[[names(rxn_rxn_networks)[i]]] = NULL
  }
  else{
    reactome_gene_metabo[[names(rxn_rxn_networks)[i]]] =
      list("gene_rxns" = pathway_genes_rxns,
           "metabo_rxns" = pathway_metabo_rxns)
  }
}

hmdb_gene_metabo = list()
for (i in 1:length(hmdb_rxn_nets)){
  pathway_name = names(hmdb_rxn_nets)[i]
  pathway_nodes = hmdb_rxn_nets[[pathway_name]] %>% unlist() %>% as.character()

  pathway_genes_rxns = hmdb_gene_pathway[[pathway_name]] %>%
    dplyr::select(PARTICIPANT, PARTICIPANT_NAME, mmu_gene) %>%
    dplyr::rename("gene" = "mmu_gene") %>%
    dplyr::mutate(rxn_id = PARTICIPANT)

  pathway_metabo_rxns = hmdb_gene_metabo_sif[[pathway_name]] %>%
    dplyr::mutate(metabo = ifelse(INTERACTION_TYPE == "consumption-controlled-by",
                                  PARTICIPANT_A,
                                  PARTICIPANT_B),
                  rxn_id = ifelse(INTERACTION_TYPE == "consumption-controlled-by",
                                  PARTICIPANT_B,
                                  PARTICIPANT_A),
                  metabo_type = ifelse(INTERACTION_TYPE == "consumption-controlled-by",
                                       "substrate",
                                       "product")) %>%
    dplyr::select(metabo, rxn_id, metabo_type)

  if(nrow(pathway_metabo_rxns) == 0 || nrow(pathway_genes_rxns) == 0){
    hmdb_gene_metabo[[pathway_name]] = NULL
  }
  else{
    hmdb_gene_metabo[[pathway_name]] =
      list("gene_rxns" = pathway_genes_rxns,
           "metabo_rxns" = pathway_metabo_rxns)
  }
}

kegg_gene_metabo = list()
for (i in 1:length(kegg_rxn_nets)){
  pathway_name = names(kegg_rxn_nets)[i]
  pathway_nodes = kegg_rxn_nets[[pathway_name]][["rxn_rxn_edges"]] %>%
    unlist() %>%
    as.character()

  rxns_info = kegg_rxn_nets[[pathway_name]][["Reactions"]] %>%
    dplyr::filter(reaction %in% pathway_nodes) %>%
    gather(key = "metabo_type", value = "kegg_cpd", -reaction, -type)

  pathway_genes_rxns = kegg_ds_genes_rxns %>%
    dplyr::filter(kegg_rxn_id %in% pathway_nodes) %>%
    dplyr::rename("rxn_id" = "kegg_rxn_id")

  pathway_metabo_rxns = kegg_ds_chebi_cpds_rxns %>%
    dplyr::filter(kegg_rxn_id %in% pathway_nodes) %>%
    dplyr::left_join(rxns_info, by = c("kegg_rxn_id" = "reaction",
                                       "kegg_cpd" = "kegg_cpd")) %>%
    dplyr::rename("metabo" = "formula") %>%
    dplyr::rename("rxn_id" = "kegg_rxn_id")

  if(nrow(pathway_metabo_rxns) == 0 || nrow(pathway_genes_rxns) == 0){
    kegg_gene_metabo[[pathway_name]] = NULL
  }
  else{
    kegg_gene_metabo[[pathway_name]] =
      list("gene_rxns" = pathway_genes_rxns,
           "metabo_rxns" = pathway_metabo_rxns)
  }
}

# Create gene-metabo pairs ---------------------------------------------------
pathway_gene_metabo = list("kegg" = kegg_gene_metabo,
                           "reactome" = reactome_gene_metabo,
                           "hmdb" = hmdb_gene_metabo)

gene_metabo_pairs_all = list()
for (i in c("kegg", "reactome", "hmdb")){
  gene_metabo_rxns_info = pathway_gene_metabo[[i]]
  pathway_gene_metabo_pairs = list()
  for (j in 1:length(gene_metabo_rxns_info)){
    pathway_name = names(gene_metabo_rxns_info)[j]

    genes = gene_metabo_rxns_info[[pathway_name]][["gene_rxns"]] %>%
      dplyr::pull(gene) %>% unique()
    metabos = gene_metabo_rxns_info[[pathway_name]][["metabo_rxns"]] %>%
      dplyr::pull(metabo) %>% unique()
    gene_metabo_pairs = expand.grid(genes, metabos)
    colnames(gene_metabo_pairs) = c("gene", "metabo")

    pathway_gene_metabo_pairs[[pathway_name]] = gene_metabo_pairs
  }
  pathway_gene_metabo_pairs = pathway_gene_metabo_pairs %>%
    dplyr::bind_rows(.id = "pathway_name")
  gene_metabo_pairs_all[[i]] = pathway_gene_metabo_pairs
}


# Calculate distance between gene-metabo pairs ----------------------------

pathway_rxn_nets = list("kegg" = lapply(kegg_rxn_nets, function(a){a[["rxn_rxn_edges"]]}),
                        "reactome" = rxn_rxn_networks,
                        "hmdb" = hmdb_rxn_nets)

pathway_rxn_nets = lapply(pathway_rxn_nets, function(a){
  a[lengths(a) != 0]
})

pathway_rxn_igraph = lapply(pathway_rxn_nets, function(x){
  lapply(x, function(y){
    graph_from_data_frame(y, directed = TRUE)
  })
})


gene_metabo_dists = list()
for (i in c("kegg", "reactome", "hmdb")){
  pairs = gene_metabo_pairs_all[[i]]
  pathways = unique(pairs$pathway_name)

  pathway_dists = list()
  for (j in pathways){
    pathway_pairs = pairs[pairs$pathway_name == j,]
    rxn_net = pathway_rxn_igraph[[i]][[j]]
    rxns_info = pathway_gene_metabo[[i]][[j]]
    pathway_dists[[j]] = Calc_gene_metabo_rxns_dist(gene_metabo_pairs = pathway_pairs,
                                                    rxns_info_list = rxns_info,
                                                    net = rxn_net)
  }
  gene_metabo_dists[[i]] = pathway_dists %>% dplyr::bind_rows()

}

all_dists = gene_metabo_dists %>% dplyr::bind_rows(.id = "pathway_source")
saveRDS(all_dists, "../../../Data/Figure 4/gene_metabo_pathway_rxn_dists.rds")

# Supp Table 18 (pathway gene-metabolite distances) is now generated by
# Analysis_scripts/Downstream/Generate_Supp_Tables.R from this saved Data/ object.



