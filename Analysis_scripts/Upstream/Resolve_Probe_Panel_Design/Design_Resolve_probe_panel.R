# Reconstructs the Resolve probe-panel gene selection (100 genes, ranked by supporting
# evidence and extended by KEGG pathway coverage). Methods-transparency only -- writes no files.

library(tidyverse)
library(readxl)

cluster_dir = Sys.getenv("SPATIAL_OMICS_CLUSTER_DIR", "")
data_dir = if (nzchar(cluster_dir)) file.path(cluster_dir, "Data/Resolve_Probe_Panel_Design") else "../../../Data/Resolve_Probe_Panel_Design"
n_additional_genes = 14

# Inputs -------------------------------------------------------------------------------------
# Per-gene HVG variance from the cell2location snRNAseq / Visium preprocessing (not reproduced
# in this repo). Each needs a `vst.variance.standardized` column and gene symbols as rownames
# (as produced by Seurat::FindVariableFeatures(..., selection.method = "vst")).
HVG_snRNAseq = readRDS(file.path(data_dir, "mouse_snRNAseq_HVGs.rds"))
HVG_visium = readRDS(file.path(data_dir, "mouse_visium_HVGs.rds"))

Resolve_probes = readxl::read_xlsx(file.path(data_dir, "Resolve_mousebrain_probe_genes.xlsx"))
Resolve_genes = Resolve_probes$GeneID

Cartana_genes_full = read.csv(file.path(data_dir, "Cartana_annotated_genes.csv"))
Cartana_genes = unique(Cartana_genes_full$Gene_symbol)

KEGG_metabolic_pathways = readRDS(file.path(data_dir, "KEGG_metabolic_pathways_clean.rds"))

# Mouse KEGG pathway -> gene membership, and Entrez -> gene symbol mapping. Live API/database
# calls, not a saved snapshot -- results may drift from what actually informed the panel design
# if KEGG or org.Mm.eg.db have been updated since (no saved version of these exists to pin to).
# keggLink() can fail with a curl HTTPS error on Windows with Bioconductor < 3.16 -- update
# Bioconductor if you hit this (https://support.bioconductor.org/p/9149679/).
mmu_KEGG_links = KEGGREST::keggLink("pathway", "mmu")
mmu_KEGG = split(names(mmu_KEGG_links), mmu_KEGG_links)

mmu_KEGG_genes = unique(unlist(mmu_KEGG))
mmu_KEGG_genes = sub(".*mmu:", "", mmu_KEGG_genes)
mapped_kegg_genes = AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = mmu_KEGG_genes,
                                           columns = "SYMBOL", keytype = "ENTREZID")

# 10x features.tsv (Ensembl ID <-> gene symbol).
symbol_ensembl = read.table(file.path(data_dir, "snRNA_Ensembl_features.tsv"))[, c(1, 2)]
colnames(symbol_ensembl) = c("Ensembl_id", "Gene_symbol")

# Manually curated first-pass metabolic-gene selection. Column 1 holds gene symbols; column 5
# flags genes excluded from that first pass.
selected_genes = readxl::read_xlsx(file.path(data_dir, "Resolve_gene_selection_draft.xlsx"),
                                    sheet = 1, col_names = FALSE)
selected_genes = selected_genes$...1[is.na(selected_genes$...5)]


# KEGG metabolic-pathway gene membership ------------------------------------------------------
entrez_to_symbol = setNames(mapped_kegg_genes$SYMBOL, mapped_kegg_genes$ENTREZID)

mmu_KEGG_metabolic_genes = mmu_KEGG[names(mmu_KEGG) %in% paste0("path:mmu", KEGG_metabolic_pathways$path_id)]
mmu_KEGG_metabolic_genes = sub(".*mmu:", "", unique(unlist(mmu_KEGG_metabolic_genes)))
mmu_KEGG_metabolic_genes = str_replace_all(mmu_KEGG_metabolic_genes, entrez_to_symbol)


# Evidence overlap across gene sources (panel-design UpSet plot) ------------------------------
venn_list = list(Resolve = Resolve_genes,
                  HVG_snRNA = rownames(HVG_snRNAseq),
                  HVG_visium = rownames(HVG_visium),
                  KEGG_metabolic = mmu_KEGG_metabolic_genes,
                  Cartana_selected = Cartana_genes)

comb_mat = ComplexHeatmap::make_comb_mat(venn_list, mode = "intersect")
ht = ComplexHeatmap::draw(ComplexHeatmap::UpSet(comb_mat, lwd = 3, pt_size = unit(5, "mm")))
comb_order = ComplexHeatmap::column_order(ht)
comb_sizes = ComplexHeatmap::comb_size(comb_mat)
ComplexHeatmap::decorate_annotation("intersection_size", {
  grid::grid.text(comb_sizes[comb_order], x = seq_along(comb_sizes),
                   y = grid::unit(comb_sizes[comb_order], "native") + grid::unit(2, "pt"),
                   default.units = "native", just = "bottom", gp = grid::gpar(fontsize = 10))
})


# Rank all candidate genes by supporting evidence ----------------------------------------------
genes_of_interest = Reduce(union, list(
  intersect(rownames(HVG_snRNAseq), mmu_KEGG_metabolic_genes),
  intersect(rownames(HVG_snRNAseq), Resolve_genes),
  intersect(rownames(HVG_visium), mmu_KEGG_metabolic_genes),
  intersect(rownames(HVG_visium), Resolve_genes),
  intersect(mmu_KEGG_metabolic_genes, Resolve_genes)
))
genes_of_interest = setdiff(genes_of_interest, Cartana_genes)

cartana_ranked = Cartana_genes_full %>%
  dplyr::select(Gene_symbol, Gene_type) %>%
  dplyr::arrange(Gene_type) %>%
  dplyr::mutate(Cartana = 1,
                HVG_snRNA = NA_real_,
                HVG_visium = NA_real_,
                KEGG_metabolic = as.numeric(Gene_type == "Metabolic"),
                Resolve = 1)

goi_ranked = purrr::map_dfr(genes_of_interest, function(gene) {
  is_metabolic = gene %in% mmu_KEGG_metabolic_genes
  tibble::tibble(
    Gene_symbol = gene,
    Gene_type = if (is_metabolic) "Metabolic" else "Cell type specific",
    HVG_snRNA = if (gene %in% rownames(HVG_snRNAseq)) HVG_snRNAseq[gene, "vst.variance.standardized"] else NA_real_,
    HVG_visium = if (gene %in% rownames(HVG_visium)) HVG_visium[gene, "vst.variance.standardized"] else NA_real_,
    KEGG_metabolic = as.numeric(is_metabolic),
    Resolve = as.numeric(gene %in% Resolve_genes),
    Cartana = 0
  )
}) %>%
  dplyr::mutate(avg_var = rowMeans(cbind(HVG_snRNA, HVG_visium), na.rm = TRUE)) %>%
  dplyr::arrange(dplyr::desc(avg_var)) %>%
  dplyr::select(-avg_var)

ranked_gene_list = dplyr::bind_rows(cartana_ranked, goi_ranked) %>%
  dplyr::left_join(symbol_ensembl, by = "Gene_symbol") %>%
  dplyr::distinct(Gene_symbol, .keep_all = TRUE) %>%
  dplyr::select(Gene_symbol, Ensembl_id, Gene_type, Cartana, HVG_snRNA, HVG_visium, KEGG_metabolic, Resolve)


# Final 100-gene metabolic panel: extend the curated selection for full KEGG pathway coverage --
kegg_pathways_for_genes = function(gene_symbols) {
  genes_entrez = mapped_kegg_genes %>%
    dplyr::filter(SYMBOL %in% gene_symbols) %>%
    dplyr::mutate(kegg_gene = paste0("mmu:", ENTREZID))

  gene_pathway_ids = purrr::map_dfr(seq_len(nrow(genes_entrez)), function(i) {
    pathway_ids = names(mmu_KEGG)[purrr::map_lgl(mmu_KEGG, ~ genes_entrez$kegg_gene[i] %in% .x)]
    tibble::tibble(kegg_gene = genes_entrez$kegg_gene[i],
                   kegg_pathway_id = sub(".*path:mmu", "", pathway_ids))
  })

  genes_entrez %>%
    dplyr::left_join(gene_pathway_ids, by = "kegg_gene") %>%
    dplyr::left_join(KEGG_metabolic_pathways, by = c("kegg_pathway_id" = "path_id")) %>%
    dplyr::filter(!is.na(pathway_names)) %>%
    dplyr::rename(Gene_symbol = SYMBOL)
}

metabolic_gene_pathways = kegg_pathways_for_genes(ranked_gene_list$Gene_symbol) %>%
  dplyr::left_join(ranked_gene_list, by = "Gene_symbol")

curated_metabolic_genes = metabolic_gene_pathways %>%
  dplyr::filter(Gene_symbol %in% selected_genes) %>%
  dplyr::pull(Gene_symbol) %>%
  unique()

# Greedily add candidate genes, in rank order, that cover at least one KEGG pathway not
# already represented by the curated set -- stops once n_additional_genes are added or
# candidates run out.
candidate_gene_order = unique(metabolic_gene_pathways$Gene_symbol)
covered_pathways = unique(metabolic_gene_pathways$pathway_names[metabolic_gene_pathways$Gene_symbol %in% curated_metabolic_genes])
additional_genes = c()

for (gene in candidate_gene_order) {
  if (length(additional_genes) == n_additional_genes) break
  if (gene %in% curated_metabolic_genes || gene %in% additional_genes) next

  gene_pathways = unique(metabolic_gene_pathways$pathway_names[metabolic_gene_pathways$Gene_symbol == gene])
  if (all(gene_pathways %in% covered_pathways)) next

  additional_genes = c(additional_genes, gene)
  covered_pathways = union(covered_pathways, gene_pathways)
}

final_100_genes = unique(c(curated_metabolic_genes, additional_genes))
final_100_panel = metabolic_gene_pathways %>%
  dplyr::filter(Gene_symbol %in% final_100_genes)
