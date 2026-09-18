
# General -----------------------------------------------------------------
map_sf_to_names <- function(df) {
  df %>%
    dplyr::left_join(metabo_molec_mapping[, c("ion", "mol_name_one_to_one")], by = c("gene" = "ion")) %>%
    dplyr::select(-gene) %>%
    dplyr::rename(gene = mol_name_one_to_one)
}


map_metabo_name_to_sf = function(metabo_name){
  metabo_molec_mapping$ion[metabo_molec_mapping$mol_name_one_to_one == metabo_name]
}
sf_to_name = function(sf){
  metabo_molec_mapping$mol_name_one_to_one[metabo_molec_mapping$ion == sf]
}

binarize_dataframe <- function(df) {
  # Apply the transformation: non-zero values become 1, zeros stay as 0
  binarized_df <- as.data.frame(lapply(df, function(x) ifelse(x != 0, 1, 0)))
  return(binarized_df)
}

subset_seurat_by_cells = function(obj, cells){
  obj@meta.data$Cell = rownames(obj@meta.data)
  obj = subset(obj, Cell %in% cells)
  return(obj)
}

jaccard_index <- function(x, y) {
  intersection <- length(intersect(x, y))
  union <- length(union(x, y))
  if (union == 0) return(NA)  # Avoid division by zero
  return(intersection / union)
}



# Regions brain -----------------------------------------------------------
rotation <- function(a){
  r <- a * pi / 180 # degrees to radians
  matrix(c(cos(r), sin(r), -sin(r), cos(r)), nrow = 2, ncol = 2)
}
get_brain_parent_depth = function(region,depth = 4){
  reg_id_path = brain_ontology$structure_id_path[brain_ontology$acronym == region]
  reg_ids = strsplit(sub(".*997/", "", reg_id_path), "/") %>% unlist()
  if (length(reg_ids) < depth){
    return(region)
  }
  else{
    parent_id = reg_ids[depth]
    parent_region = brain_ontology$acronym[brain_ontology$id == parent_id]
    return(parent_region)
  }
}
get_children_regions = function(ontology_df, ROI_acronym){
  id = ontology_df$id[ontology_df$acronym == ROI_acronym]
  child_paths = ontology_df$structure_id_path[str_detect(ontology_df$structure_id_path,
                                                         paste0(id, "/"))]
  children = sub(paste0(".*", id), "", child_paths)
  children = strsplit(children, "/") %>%
    unlist() %>%
    unique()
  children = children[children != ""] %>% as.integer()
  children_regions = ontology_df$acronym[ontology_df$id %in% children]
  return(children_regions)
}


# MISTY -------------------------------------------------------------------
Get_interaction_importance = function(misty.results, view, cutoff = 1,
                                      trim = -Inf, trim.measure = c(
                                        "gain.R2", "multi.R2", "intra.R2",
                                        "gain.RMSE", "multi.RMSE", "intra.RMSE"
                                      ),clean = FALSE){
  trim.measure.type <- match.arg(trim.measure)

  assertthat::assert_that(("importances.aggregated" %in% names(misty.results)),
                          msg = "The provided result list is malformed. Consider using collect_results()."
  )

  assertthat::assert_that(("improvements.stats" %in% names(misty.results)),
                          msg = "The provided result list is malformed. Consider using collect_results()."
  )

  assertthat::assert_that((view %in%
                             (misty.results$importances.aggregated %>% dplyr::pull(view))),
                          msg = "The selected view cannot be found in the results table."
  )

  inv <- sign((stringr::str_detect(trim.measure.type, "gain") |
                 stringr::str_detect(trim.measure.type, "RMSE", negate = TRUE)) - 0.5)

  targets <- misty.results$improvements.stats %>%
    dplyr::filter(
      measure == trim.measure.type,
      inv * mean >= inv * trim
    ) %>%
    dplyr::pull(target)


  plot.data <- misty.results$importances.aggregated %>%
    dplyr::filter(view == !!view,
                  Target %in% targets,
                  Importance >= cutoff)

  if (clean) {
    clean.predictors <- plot.data %>%
      dplyr::mutate(Importance = Importance * (Importance >= cutoff)) %>%
      dplyr::group_by(Predictor) %>%
      dplyr::summarize(total = sum(Importance, na.rm = TRUE)) %>%
      dplyr::filter(total > 0) %>%
      dplyr::pull(Predictor)
    clean.targets <- plot.data %>%
      dplyr::mutate(Importance = Importance * (Importance >= cutoff)) %>%
      dplyr::group_by(Target) %>%
      dplyr::summarize(total = sum(Importance, na.rm = TRUE)) %>%
      dplyr::filter(total > 0) %>%
      dplyr::pull(Target)
    plot.data.clean <- plot.data %>%
      dplyr::filter(
        Predictor %in% clean.predictors,
        Target %in% clean.targets
      )
  } else {
    plot.data.clean <- plot.data
  }
  return(plot.data.clean)
}

prepare_views_misty = function(predictor_data, target_data, pred_view_name){
  common_cells = intersect(rownames(predictor_data),
                           rownames(target_data))

  target_view = mistyR::create_initial_view(target_data[common_cells,])

  predictor_view = mistyR::create_view(pred_view_name,
                               predictor_data[common_cells,],
                               pred_view_name)

  combined_views = target_view %>%
    mistyR::add_views(predictor_view)

  return(combined_views)

}

create_reg_cell_type_composition = function(ct_comp, reg_comp,
                                            intersect_data_long, min_cells = 100){
  filt = intersect_data_long %>%
    dplyr::filter(n_cells >= min_cells)

  ct_region_cells = intersect(rownames(reg_comp),
                              rownames(ct_comp))

  combined = list()
  for (i in 1:nrow(filt)){
    reg_name = filt$Region[i]
    ct_name = filt$Cell_type[i]

    res = data.frame(v1 = ct_region_cells,
                     v2 = ct_comp[ct_region_cells,ct_name] * reg_comp[ct_region_cells, reg_name])
    colnames(res) = c("Cell", paste0(ct_name, ".", reg_name))

    combined[[i]] = res %>%
      column_to_rownames("Cell") %>%
      as.data.frame()
  }

  combined = dplyr::bind_cols(combined)

  combined = combined[rowSums(combined) > 0,]

  return(combined)
}

# Metabo-gene correlation -------------------------------------------------
extract_cor_data = function(seurat_obj, data_type,
                            rna_feats = NULL, metabo_feats = NULL,
                            metabo_LOD = 100){

  metabo_key = Key(seurat_obj[["Metabo"]])
  rna_key = Key(seurat_obj[["RNA"]])

  if(!is.null(rna_feats)){
    rna_feats = paste0(rna_key, rna_feats)
  }
  else{
    rna_feats = paste0(rna_key, rownames(seurat_obj[["RNA"]]))
  }

  if(!is.null(metabo_feats)){
    metabo_feats = paste0(metabo_key, metabo_feats)
  }
  else{
    metabo_feats = paste0(metabo_key, rownames(seurat_obj[["Metabo"]]))
  }

  cor_data = FetchData(seurat_obj, vars = c(rna_feats, metabo_feats),
                       slot = data_type)
  counts_data = FetchData(seurat_obj, vars = c(rna_feats, metabo_feats),
                          slot = "counts")

  cor_data[,metabo_feats][counts_data[,metabo_feats] < metabo_LOD] = 0
  if (data_type == "counts"){
    cor_data = cor_data %>%
      mutate_all(as.integer)
  }

  return(cor_data)

}

contextualize_seurat_obj <- function(seurat_obj,
                                     metadata_filters,
                                     min_cells = 30) {
  # Check that all metadata columns exist in the Seurat object
  seurat_obj@graphs = list()

  if(length(metadata_filters) == 0){
    return(seurat_obj)
  }

  missing_cols <- setdiff(names(metadata_filters), colnames(seurat_obj@meta.data))
  if (length(missing_cols) > 0) {
    stop(paste("The following metadata columns are missing:", paste(missing_cols, collapse = ", ")))
  }

  # Start with all cell names
  selected_cells <- colnames(seurat_obj)

  # Iteratively filter cell names based on each metadata condition
  for (col_name in names(metadata_filters)) {
    value <- metadata_filters[[col_name]]
    selected_cells <- intersect(
      selected_cells,
      rownames(seurat_obj@meta.data[seurat_obj@meta.data[[col_name]] == value, ])
    )
  }

  # Check if number of selected cells meets the minimum threshold
  if (length(selected_cells) >= min_cells) {
    return(subset(seurat_obj, cells = selected_cells))
  } else {
    message(paste("Only", length(selected_cells), "cells match the criteria. Minimum required is", min_cells, ". Returning NULL."))
    return(NULL)
  }
}

custom_if_dropout_scimpute = function(mat, ncores, dthre){
  mat[mat <= log10(1.01)] = log10(1.01)
  pa = scLink:::get_mix_parameters(t(mat), ncores = ncores)
  I = ncol(mat)
  J= nrow(mat)
  ### calculate cell-wise dropout rate
  droprate = sapply(1:I, function(i) {
    if(is.na(pa[i,1])) return(rep(1,J)) #Changed from original rep(0,J) to rep(1,J)
    wt = scLink:::calculate_weight(mat[, i], pa[i, ])
    return(wt[, 1])
  })
  dropind = 1* (droprate > dthre)
  return(dropind)
}

sclink_est_scimpute_custom = function(count, cor.p,
                                      corr_method = "pearson",
                                      thre_no = 20,
                                      dthre = 0.9, ncores){
  count_dr = count
  x = log10(10^count - 1 + 1.01)
  dropind = custom_if_dropout_scimpute(x, ncores = ncores, dthre)
  count_dr[dropind == 1] = NA
  count_NA_inds = count_dr
  cor_complete = cor(count_dr,method = corr_method,
                     use = "pairwise.complete.obs")
  cor_complete[is.na(cor_complete)] = cor.p[is.na(cor_complete)]

  count_dr = 1 * (!is.na(count_dr))
  nobs = t(count_dr) %*% count_dr
  cor_complete[nobs < thre_no] = cor.p[nobs < thre_no]
  return(list("cor_res" = cor_complete,
              "cor_input" = count_NA_inds))
}

scLink_n_obs_corr_pivot_long = function(impute_ind_mat){
  count_dr = 1 * (!is.na(impute_ind_mat))
  nobs = t(count_dr) %*% count_dr %>%
    as.data.frame() %>%
    rownames_to_column("feat_A") %>%
    pivot_longer(cols = -feat_A,
                 names_to = "feat_B",
                 values_to = "n_cells_cor") %>%
    dplyr::filter(!is.na(n_cells_cor),
                  feat_A != feat_B)
  return(nobs)
}

sclink_cor_custom = function(expr, ncores, corr_method = "pearson",
                             nthre = 20, dthre = 0.9){
  cor_pearson = cor(expr,method = corr_method)
  x.q = apply(expr, 2, sd)
  cor_pearson.c = sclink_est_scimpute_custom(expr, ncores = ncores,
                                             cor.p = cor_pearson,
                                             thre_no = nthre, dthre = dthre)
  return(cor_pearson.c)

}

calc_rna_metabo_corr = function(seurat_obj,
                                corr_method,
                                metadata_filters,
                                rna_feats = NULL, metabo_feats = NULL,
                                min_nonzero_cells = 20,
                                feat_nonzero_prop = 0.9,
                                metabo_LOD = 100,
                                nonzero_only = F){

  context_seurat = contextualize_seurat_obj(seurat_obj = seurat_obj,
                                            metadata_filters = metadata_filters,
                                            min_cells = min_nonzero_cells)
  if(is.null(context_seurat)){
    return(NULL)
  }

  cor_data = extract_cor_data(seurat_obj = context_seurat, data_type = "data",
                              rna_feats = rna_feats,
                              metabo_feats = metabo_feats,
                              metabo_LOD = metabo_LOD) %>% as.matrix()

  if (nonzero_only){
    cor_data[cor_data == 0] = NA
  }
  scLink_cor_res = sclink_cor_custom(expr = cor_data,ncores = 1,
                                     corr_method = corr_method,
                                     nthre = min_nonzero_cells,
                                     dthre = feat_nonzero_prop)

  scLink_cor_res$cor_res[is.na(scLink_cor_res$cor_res)] = 0

  cor_input = scLink_cor_res$cor_input
  n_cells_obs = scLink_n_obs_corr_pivot_long(cor_input)

  scLink_cor_df = scLink_cor_res$cor_res %>%
    as.data.frame() %>%
    rownames_to_column("feat_A") %>%
    pivot_longer(cols = -feat_A,
                 names_to = "feat_B",
                 values_to = "coef") %>%
    dplyr::filter(!is.na(coef),
                  feat_A != feat_B) %>%
    dplyr::left_join(n_cells_obs)

  scLink_cor_df$feat_A_type = ifelse(stringr::str_starts(scLink_cor_df$feat_A, "metabo"),
                                     "Metabo",
                                     "RNA")
  scLink_cor_df$feat_B_type = ifelse(stringr::str_starts(scLink_cor_df$feat_B, "metabo"),
                                     "Metabo",
                                     "RNA")

  return(list("cor_df" = scLink_cor_df,
              "cor_input" = cor_input))
}

get_cor_cells_input = function(impute_ind_mat,rna_feat, metabo_feat){

  count_dr = 1 * (!is.na(impute_ind_mat))
  cells = count_dr[,c(rna_feat, metabo_feat)] %>% rowSums()
  cells = cells[cells == 2] %>% names()

  return(cells)
}


# MCP DIALOGUE -------------------------------------------------
create_celltype_dialogue = function(seurat_obj,
                                    niches_res,
                                    assay = "Metabo",
                                    cell_type){
  cell_type_name = cell_type

  tpm = seurat_obj@assays[[assay]]$data
  counts = seurat_obj@assays[[assay]]$counts

  ct_cells = seurat_obj@meta.data %>%
    dplyr::filter(Cell_type_coarse == cell_type_name) %>%
    rownames()

  ct_niches = niches_res[niches_res$Cell %in% ct_cells,]

  tpm = tpm[,ct_niches$Cell] %>% as.data.frame()
  counts = counts[,ct_niches$Cell] %>% as.data.frame()

  samples = ct_niches$Niche

  ct_meta = metabo_RNA_seurat@images$image@coordinates %>%
    dplyr::select(x,y) %>%
    dplyr::mutate(x = -1 * x,
                  y = -1 * y) %>%
    rownames_to_column("Cell") %>%
    dplyr::distinct() %>%
    dplyr::filter(Cell %in% ct_niches$Cell) %>%
    dplyr::left_join(ct_niches) %>%
    dplyr::rename("coor.Y" = "x",
                  "coor.X" = "y",
                  ) %>%
    column_to_rownames("Cell")

  ct_meta = ct_meta[ct_niches$Cell,]

  n_count = colSums(counts) %>%
    log() %>% scale() %>% as.vector()

  dialoge_ct = DIALOGUE::make.cell.type(name = cell_type_name,
                              tpm = tpm,
                              samples = samples,
                              X = t(tpm),
                              metadata = ct_meta,
                              cellQ = n_count)
  return(dialoge_ct)
}

calculate_environment_score = function(ra_res){
  n_cell_types = length(ra_res) - 2
  all_samples = lapply(ra_res[1:n_cell_types], function(x){x@samples}) %>%
    unlist() %>%
    unique()

  res_per_sample = list()
  for (i in all_samples){
    ct_cells_sample = lapply(ra_res[1:n_cell_types], function(x){
      x@cells[x@samples == i]
    })
    ct_cells_sample = ct_cells_sample[lengths(ct_cells_sample) != 0]
    if (length(ct_cells_sample) < 2){
      next()
    }
    ct_res = list()
    for (ct in 1:length(ct_cells_sample)){
      ct_name = names(ct_cells_sample)[ct]
      ct_scores = ra_res[[ct_name]]@scores[ct_cells_sample[[ct]],]
      other_cell_types = names(ct_cells_sample)[-ct]
      other_ct_scores = lapply(ra_res[other_cell_types], function(x){
        n_cells = length(x@cells[x@samples == i])
        other_score = x@scores[x@cells[x@samples == i],]
        if(n_cells == 1){
          other_score = other_score %>% data.frame() %>% t()
        }
        colSums(other_score) / n_cells
      })

      expected = other_ct_scores %>%
        dplyr::bind_rows() %>%
        colSums()

      if(length(ct_cells_sample[[ct]]) == 1){
        env_score = sweep(t(ct_scores), 2, expected, FUN = "-") %>% abs()
      }
      else{
        env_score = sweep(ct_scores, 2, expected, FUN = "-") %>% abs()
      }

      Es_ij = rowSums(env_score) * -1
      Es_ij = data.frame(Es_ij) %>%
        rownames_to_column("Cell") %>%
        dplyr::mutate(sample = i) %>%
        dplyr::mutate(cell_type = ct_name) %>%
        dplyr::relocate(sample,cell_type, .after = Cell)

      if(length(ct_cells_sample[[ct]]) == 1){
        Es_ij$Cell = ct_cells_sample[[ct]] %>% as.character()
      }

      ct_res[[ct_name]] = Es_ij
    }
    res_per_sample[[i]] = ct_res %>% dplyr::bind_rows()
  }
  final = res_per_sample %>% dplyr::bind_rows()
  return(final)
}

Merge_MCP_signature = function(sig_per_ct){
  merged_res = list()
  counter = 0
  for (ct in 1:length(sig_per_ct)){
    for (i in 1:length(sig_per_ct[[ct]])){
      counter = counter + 1
      MCP_direction = strsplit(names(sig_per_ct[[ct]])[i], "[.]") %>% unlist()
      merged_res[[counter]] = data.frame(Cell_type = names(sig_per_ct)[ct],
                                         Metabos = sig_per_ct[[ct]][[i]],
                                         MCP = MCP_direction[1],
                                         Direction = MCP_direction[2])
    }
  }
  merged_res = merged_res %>% dplyr::bind_rows()
  merged_res$Direction = ifelse(merged_res$Direction == "up", 1, -1)
  return(merged_res)
}

custom_D3_DIALOGUE = function(rA_list,D2_result){
  R = D2_result
  full_rA = rA_list
  n_cell_types = length(rA_list)

  R$gene.pval<-lapply(R$cell.types, function(x){DIALOGUE:::DLG.multi.get.gene.pval(x,R)})
  names(R$gene.pval)<-R$cell.types

  full_rA<-lapply(full_rA,function(r){
    r<-DIALOGUE:::DLG.find.scoring(r,R)
    return(r)
  })

  full_rA$scores <-lapply(full_rA,function(r1){
    X<-cbind.data.frame(r1@scores,samples = r1@samples,
                        cells = r1@cells, cell.type = r1@name,
                        r1@metadata)
    return(X)})

  full_rA$scores_samples <-lapply(full_rA[1:n_cell_types],function(r1){
    X<-r1@scoresAv
    rownames(X) = paste0(r1@name, ";", rownames(X))
    X = as.data.frame(X)
    return(X)})

  cell.types = names(full_rA)[1:n_cell_types]
  R$pref<-list()
  idx<-unique(get.strsplit(names(R$sig[[1]]),".",1))
  pairs1<-t(combn(cell.types,2))

  for(i in 1:nrow(pairs1)){
    x1<-pairs1[i,1];x2<-pairs1[i,2]
    x<-paste0(x1,".vs.",x2)
    r1<-full_rA[[x1]];r2<-full_rA[[x2]]
    r<-DIALOGUE:::DLG.get.OE(r1,r2,plot.flag = F,compute.scores = F)
    r1<-r$r1;r2<-r$r2
    idx<-intersect(DIALOGUE:::get.abundant(r1@samples),DIALOGUE:::get.abundant(r2@samples))
    R$pref[[x]]<-cbind.data.frame(R = diag(cor(r1@scoresAv[idx,],r2@scoresAv[idx,])),
                                  hlm = DIALOGUE:::DLG.hlm.pval(r1,r2,formula = R$frm))
  }
  R$gene.pval<-lapply(full_rA[1:n_cell_types],function(r1) r1@gene.pval)
  R$sig1<-lapply(full_rA[1:n_cell_types],function(r1) r1@sig$sig1)
  R$sig2<-lapply(full_rA[1:n_cell_types],function(r1) r1@sig$sig2)

  R$scores <-lapply(full_rA[1:n_cell_types],function(r1){
    X<-cbind.data.frame(r1@scores,samples = r1@samples,
                        cells = r1@cells, cell.type = r1@name,
                        r1@metadata)
    return(X)})

  R$cca.fit<-plyr::laply(R$cell.types,function(x) diag(cor(R$cca.scores[[x]],R$scores[[x]][,1:R$k["DIALOGUE"]])))
  rownames(R$cca.fit)<-R$cell.types

  return(list("R" = R,
              "full_RA_res" = full_rA))


}

MCP_signature_to_df = function(sig_list, min_celltype = 4){

  df = lapply(sig_list, function(a){
    enframe(a, name = "MCP_reg", value = "metabolite") %>%
      unnest(metabolite)}) %>%
    dplyr::bind_rows(.id = "cell_type") %>%
    tidyr::separate(MCP_reg, into = c("MCP", "Direction"),
                    sep = "[.]")

  df_n_ct = df %>%
    dplyr::group_by(metabolite, MCP, Direction) %>%
    dplyr::summarise(n_CT = n()) %>%
    dplyr::filter(n_CT >= min_celltype)

  final = df %>%
    dplyr::left_join(df_n_ct) %>%
    dplyr::filter(!is.na(n_CT))

  return(final)
}
# Pathway reaction networks -----------------------------------------------
Get_pathways_of_rxns <- function(rxn_id) {
  base_url <- "https://reactome.org/ContentService/data/entity/"
  full_url <- paste0(base_url, rxn_id, "/componentOf")

  tryCatch({
    response <- httr::GET(full_url)

    if (httr::status_code(response) != 200) {
      return(NULL)
    }

    result <- httr::content(response, as = "parsed", type = "application/json")
    pathway_info = data.frame(pathway_id = unlist(result[[1]]$stIds),
                              pathway_name = unlist(result[[1]]$names),
                              rxn_id = rxn_id)
    return(pathway_info)

  }, error = function(e) {
    return(NULL)
  })
}
Get_rxns_in_pathway = function(pathway_id){
  base_url <- "https://reactome.org/ContentService/data/pathway/"
  full_url <- paste0(base_url, pathway_id, "/containedEvents")

  tryCatch({
    response <- httr::GET(full_url)

    if (httr::status_code(response) != 200) {
      return(NULL)
    }

    result <- httr::content(response, as = "parsed", type = "application/json")
    rxns = lapply(result, function(x){
      ifelse(x$schemaClass == "Reaction",x$stId, NA)
    }) %>% unlist()
    rxns = rxns[!is.na(rxns)]
    final = data.frame(rxn_id = rxns, pathway_id = pathway_id)
    return(final)

  }, error = function(e) {
    return(NULL)
  })
}
Get_rxns_preceding = function(rxn_id){
  base_url <- "https://reactome.org/ContentService/data/query/enhanced/"
  full_url <- paste0(base_url, rxn_id)

  tryCatch({
    response <- httr::GET(full_url)

    if (httr::status_code(response) != 200) {
      return(NULL)
    }

    result <- httr::content(response, as = "parsed", type = "application/json")

    from_rxns = lapply(result$precedingEvent, function(x){
      x$stId}) %>% unlist()

    final = data.frame(from = from_rxns, to = rxn_id)

    return(final)

  }, error = function(e) {
    return(NULL)
  })
}

kegg_add_reversible_edges <- function(edges_df, reactions_df) {
  # Merge to get type info for 'from' and 'to' reactions
  edges_with_types <- merge(edges_df, reactions_df, by.x = "from", by.y = "reaction")
  colnames(edges_with_types)[colnames(edges_with_types) == "type"] <- "from_type"

  edges_with_types <- merge(edges_with_types, reactions_df, by.x = "to", by.y = "reaction")
  colnames(edges_with_types)[colnames(edges_with_types) == "type"] <- "to_type"

  # Find rows where both from and to are reversible
  reversible_edges <- edges_with_types[
    edges_with_types$from_type == "reversible" &
      edges_with_types$to_type == "reversible", c("from", "to")
  ]

  # Create reversed edges
  reversed_edges <- reversible_edges
  reversed_edges$from <- reversible_edges$to
  reversed_edges$to <- reversible_edges$from

  # Combine original and reversed edges
  updated_edges <- rbind(edges_df, reversed_edges)

  return(updated_edges)
}
kegg_pathway2rxn_network = function(pathway_id){
  url = paste0("https://rest.kegg.jp/get/", pathway_id, "/kgml")
  status_code = httr::GET(url) %>% httr::status_code()

  if(status_code != 200){
    return(NULL)
  }

  tmp = tempfile()
  utils::download.file(
    url,
    destfile = tmp)

  parsed_kgml = KEGGgraph::parseKGML(tmp)
  pathway_graph = KEGGgraph::KEGGpathway2Graph(parsed_kgml,genesOnly = F)
  edge_data = KEGGgraph::getKEGGedgeData(pathway_graph)
  rxns = KEGGgraph::getReactions(parsed_kgml)

  if(length(rxns) == 0){
    return(NULL)
  }

  if(length(edge_data) == 0){
    return(NULL)
  }


  rxns_info = lapply(rxns, function(x){
    data.frame(reaction = KEGGgraph::getName(x),
               type = KEGGgraph::getType(x),
               substrate = KEGGgraph::getSubstrate(x) %>% paste(collapse = ";"),
               product = KEGGgraph::getProduct(x) %>% paste(collapse = ";"))}) %>%
    dplyr::bind_rows()

  rxn_type_info = rxns_info %>%
    dplyr::select(reaction, type) %>%
    dplyr::distinct()

  rxn_rxn_edges = lapply(edge_data, function(x){
    KEGGgraph::getEntryID(x) %>%
      t() %>%
      data.frame()}) %>%
    dplyr::bind_rows() %>%
    dplyr::filter(Entry1ID %in% rxns_info$reaction,
                  Entry2ID %in% rxns_info$reaction)
  colnames(rxn_rxn_edges) = c("from", "to")

  final_edges = kegg_add_reversible_edges(edges_df = rxn_rxn_edges,
                                          reactions_df = rxn_type_info)

  return(list("Reactions" = rxns_info,
              "rxn_rxn_edges" = final_edges))
}

#This one is for RHEA
calc_shortest_path_and_distance = function(net, root_node,
                                           metabolites_of_interest = NULL){
  enzyme_reaction_node = root_node
  network = net


  reachable_nodes <- bfs(network,
                         root = enzyme_reaction_node, mode = "all")$order
  reachable_nodes <- reachable_nodes[!is.na(reachable_nodes)]
  subnetwork <- induced_subgraph(network, vids = reachable_nodes)

  metabolites <- V(subnetwork)$name[V(subnetwork)$type == "metabolite"]

  if(!is.null(metabolites_of_interest)){
    metabolites = intersect(metabolites, metabolites_of_interest)
  }

  shortest_paths <- list()
  distances <- c()

  for (i in c("upstream", "downstream")) {
    for (metabolite in metabolites){
      if (i == "upstream"){
        path <- shortest_paths(network, from = metabolite, to = enzyme_reaction_node,
                               output = "both", mode = "out")
        distance <- distances(network, v = metabolite, to = enzyme_reaction_node)

        if(length(path$epath[[1]]) > 0){
          shortest_paths[[i]][[metabolite]] <- path
          distances[[i]][metabolite] <- distance
        }
      }
      else{
        path <- shortest_paths(network, from = enzyme_reaction_node, to = metabolite,
                               output = "both", mode = "out")
        distance <- distances(network, v = enzyme_reaction_node, to = metabolite)

        if(length(path$epath[[1]]) > 0){
          shortest_paths[[i]][[metabolite]] <- path
          distances[[i]][metabolite] <- distance
        }
      }
    }
  }

  return(list("paths" = shortest_paths,
              "distances" = distances))
}



# Metabo-gene distance reaction networks ----------------------------------
get_path_length <- function(net, source_node, target_node) {
  # Build directed graph from edge list

  # Check if a path exists
  if (suppressWarnings(!are.connected(net, source_node, target_node))) {
    return(NA)  # No path exists
  }
  else{
    path_length = 1
  }

  # Get the shortest path length (number of edges)
  # path_length <- lengths(shortest_paths(net,
  #                                       from = source_node,
  #                                       to = target_node,
  #                                       output = "epath")$epath[[1]])

  # Return path length (0 if same node, >0 if a path exists)
  return(path_length)
}
Calc_gene_metabo_rxns_dist = function(gene_metabo_pairs,
                                      rxns_info_list,
                                      net){
  if(!"metabo_type" %in% colnames(rxns_info_list[["metabo_rxns"]])){
    rxns_info_list[["metabo_rxns"]] = rxns_info_list[["metabo_rxns"]] %>%
      dplyr::mutate(metabo_type = NA)
  }
  rxns_info_list[["metabo_rxns"]]  = rxns_info_list[["metabo_rxns"]] %>%
    dplyr::select(metabo, rxn_id, metabo_type) %>%
    dplyr::rename("metabo_rxn_id" = "rxn_id")

  rxns_info_list[["gene_rxns"]]  = rxns_info_list[["gene_rxns"]] %>%
    dplyr::select(gene, rxn_id) %>%
    dplyr::rename("gene_rxn_id" = "rxn_id")

  mols_rxns = gene_metabo_pairs %>%
    dplyr::left_join(rxns_info_list[["gene_rxns"]], by = "gene") %>%
    dplyr::left_join(rxns_info_list[["metabo_rxns"]], by = "metabo") %>%
    dplyr::distinct()

  path_lengths = integer(nrow(mols_rxns))
  for (i in 1:nrow(mols_rxns)){
    if(!is.na(mols_rxns$metabo_type[i])){
      if(mols_rxns$metabo_type[i] == "substrate"){
        source = mols_rxns$metabo_rxn_id[i]
        target = mols_rxns$gene_rxn_id[i]
      }
      else{
        source = mols_rxns$gene_rxn_id[i]
        target = mols_rxns$metabo_rxn_id[i]
      }
      if(!source %in% V(net)$name || !target %in% V(net)$name){
        path_lengths[i] = NA
      }
      else{
        if(source == target){
          path_lengths[i] = 0
        }
        else{
          dist = get_path_length(net = net,
                                 source_node = source,
                                 target_node = target)
          path_lengths[i] = dist
        }
      }
    }
    else{
      up_source = mols_rxns$metabo_rxn_id[i]
      up_target = mols_rxns$gene_rxn_id[i]

      down_source = mols_rxns$gene_rxn_id[i]
      down_target = mols_rxns$metabo_rxn_id[i]

      if(!up_source %in% V(net)$name || !up_target %in% V(net)$name){
        path_lengths[i] = NA
      }
      else{
        if(up_source == up_target){
          path_lengths[i] = 0
        }
        else{
          up_dist = get_path_length(net = net,
                                    source_node = up_source,
                                    target_node = up_target)
          down_dist = get_path_length(net = net,
                                      source_node = down_source,
                                      target_node = down_target)
          dist = min(up_dist, down_dist, na.rm = T)
          if(is.infinite(dist)){
            dist = NA
          }
          path_lengths[i] = dist
        }
      }
    }
  }

  mols_rxns$dist = path_lengths
  return(mols_rxns)

}

# Compare DE to DIALOGUE -----------------------------------------------------
Get_metabotype_MCP_markers = function(dialogue_results,
                                      MCP_metabotype_cells,
                                      metabotype_DE,
                                      metabotype,
                                      include_DE_overlap = F,
                                      min_MCP_metabotype_prop = 0.8){
  MCPs_to_consider = conf_mat_MCP_metabotype %>%
    dplyr::filter(Metabotype == metabotype,
                  prop_common >= min_MCP_metabotype_prop) %>%
    dplyr::pull(MCP)
  if(length(MCPs_to_consider) == 0){
    return(NULL)
  }
  else{
    all_up_markers = c()
    MCP_metabotyp_common_cells = c()
    for (mcp in MCPs_to_consider){
      if (! mcp %in% dialogue_sig$MCP){
        return(NULL)
      }
      sig_up = dialogue_sig %>%
        dplyr::filter(MCP == mcp,
                      n_CT == 5) %>%
        dplyr::pull(metabolite)
      MCP_cells = MCP_metabotype_cells %>%
        dplyr::filter(MCP == mcp,
                      Metabotype == metabotype) %>%
        dplyr::pull(Cell)
      all_up_markers = c(all_up_markers, sig_up)
      MCP_metabotyp_common_cells = c(MCP_metabotyp_common_cells,
                                     MCP_cells)
    }

    MCP_metabotyp_common_cells = unique(MCP_metabotyp_common_cells)

    DE_metabotype_markers = metabotype_DE %>%
      dplyr::filter(avg_log2FC >= 1,
                    p_val_adj < 0.05,
                    cluster == metabotype) %>%
      dplyr::pull(gene)

    if(!include_DE_overlap){
      final_markers = setdiff(all_up_markers, DE_metabotype_markers)
    }
    else{
      final_markers = unique(all_up_markers)
    }

    only_metabotype_markers = setdiff(DE_metabotype_markers,all_up_markers)

    MCP_metabotype_expr = list()
    only_metabotype_expr = list()
    main_celltypes = c("Astrocytes", "Excit", "Inhib",
                       "Oligodendrocytes", "Microglia")
    for (m in final_markers){
      marker_expr = list()
      for (i in main_celltypes){
        common_cells = intersect(MCP_metabotyp_common_cells,
                                 colnames(dialogue_results[[i]]@tpm))

        if(length(common_cells) == 0){
          next()
        }

        expr_vals = dialogue_results[[i]]@tpm[m,common_cells]
        if(median(as.numeric(expr_vals)) < 1){
          next()
        }
        marker_expr[[i]] = data.frame(Cell_type = i,metabolite = m,
                                exp = as.numeric(expr_vals))
      }
      if (length(marker_expr) == 0){
        next()
      }
      marker_expr = marker_expr %>% dplyr::bind_rows()

      MCP_metabotype_expr[[m]] = marker_expr
    }
    MCP_metabotype_expr = MCP_metabotype_expr %>% dplyr::bind_rows()

    for (m in only_metabotype_markers){
      marker_expr = list()
      for (i in main_celltypes){
        common_cells = intersect(MCP_metabotyp_common_cells,
                                 colnames(dialogue_results[[i]]@tpm))

        if(length(common_cells) == 0){
          next()
        }

        expr_vals = dialogue_results[[i]]@tpm[m,common_cells]
        marker_expr[[i]] = data.frame(Cell_type = i,metabolite = m,
                                exp = as.numeric(expr_vals))
      }
      if (length(marker_expr) == 0){
        next()
      }
      marker_expr = marker_expr %>% dplyr::bind_rows()

      only_metabotype_expr[[m]] = marker_expr
    }
    MCP_metabotype_expr = MCP_metabotype_expr %>% dplyr::bind_rows()
    only_metabotype_expr = only_metabotype_expr %>% dplyr::bind_rows()

    n_df <- MCP_metabotype_expr %>%
      dplyr::group_by(Cell_type, metabolite) %>%
      dplyr::summarise(n = n(), .groups = "drop") %>%
      dplyr::select(Cell_type,n) %>%
      dplyr::distinct()

    n_df_metabotype = only_metabotype_expr %>%
      dplyr::group_by(Cell_type, metabolite) %>%
      dplyr::summarise(n = n(), .groups = "drop") %>%
      dplyr::select(Cell_type,n) %>%
      dplyr::distinct()

    MCP_metabotype_expr = MCP_metabotype_expr %>%
      dplyr::left_join(n_df) %>%
      dplyr::mutate(Cell_type = paste0(Cell_type, "( ", n, ")"))

    only_metabotype_expr = only_metabotype_expr %>%
      dplyr::left_join(n_df_metabotype) %>%
      dplyr::mutate(Cell_type = paste0(Cell_type, "( ", n, ")"))

    sel_markers = unique(MCP_metabotype_expr$metabolite)

    p1 = ggboxplot(MCP_metabotype_expr,
                   x = "metabolite",
                   y = "exp",
                   fill = "Cell_type",width = 0.5,
                   outlier.shape = NA) +
      rotate_x_text(angle = 45)

    metabotype_plot_seurat = metabo_RNA_seurat
    Idents(metabotype_plot_seurat) = "Metabotype"

    metabotype_cells = rownames(metabo_RNA_seurat@meta.data)[metabo_RNA_seurat@meta.data$Metabotype == metabotype]


    sel_cells = list("Met.2" = metabotype_cells,
                     "MCP" = MCP_metabotyp_common_cells)
    names(sel_cells) = c(metabotype, paste(MCPs_to_consider, collapse = ";"))

    p2 = SpatialDimPlot(object = metabotype_plot_seurat,
                        cells.highlight = sel_cells,
                        cols.highlight = c("darkred",
                                           "grey"),
                        stroke = NA,facet.highlight = T,
                        pt.size.factor = 2)

    p3 = ggboxplot(only_metabotype_expr,
                   x = "metabolite",
                   y = "exp",
                   fill = "Cell_type",width = 0.7,
                   outlier.shape = NA) +
      rotate_x_text(angle = 45)

    return(list("markers" = sel_markers,
                "Expr_data" = MCP_metabotype_expr,
                "Expr_data_metabotype_only" = only_metabotype_expr,
                "plots" = list("expr_boxplot" = p1,
                               "spatial" = p2,
                               "expr_boxplot_metabotype_only" = p3)))

  }


}



# Compare DE to MISTY -----------------------------------------------------

compare_DE_misty = function(DE_res, misty_res,
                            misty_ion_name_mapping,
                            LFC_cutoff = 1,
                            padj_cutoff = 0.05,
                            min_pct = 0.5,
                            importance_cutoff = 2,
                            trim.measure = "gain.RMSE",
                            trim = 10){
  DE_col_names = c(Metabolite = "gene", Cluster = "cluster",
                   Score = "avg_log2FC")
  Misty_col_names = c(Metabolite = "Predictor", Cluster = "Target",
                      Score = "Importance")

  direction_data = DE_res %>%
    dplyr::select(gene,cluster,avg_log2FC) %>%
    dplyr::rename(all_of(DE_col_names)) %>%
    dplyr::mutate(Direction = ifelse(sign(Score) == 1, "Up", "Down")) %>%
    dplyr::select(-Score) %>%
    dplyr::distinct() %>%
    dplyr::mutate(Metabolite = gsub("[+-]", "_", Metabolite))

  DE_res_scores = DE_res %>%
    dplyr::filter(avg_log2FC >= LFC_cutoff,
                  p_val_adj < padj_cutoff,
                  pct.1 >= min_pct) %>%
    dplyr::select(gene,cluster,avg_log2FC) %>%
    dplyr::rename(all_of(DE_col_names)) %>%
    dplyr::select(-Score)

  misty_res_scores = Get_interaction_importance(misty.results = misty_res,
                                                view = "Metabo",
                                                cutoff = importance_cutoff,
                                                trim.measure = trim.measure,
                                                trim = trim, clean = T) %>%
    dplyr::left_join(misty_ion_name_mapping, by = "Predictor") %>%
    dplyr::select(-Predictor) %>%
    dplyr::rename("Predictor" = "ion") %>%
    dplyr::select(Predictor, Target, Importance) %>%
    dplyr::rename(all_of(Misty_col_names)) %>%
    dplyr::left_join(direction_data) %>%
    # dplyr::filter(Direction == "Up") %>%
    dplyr::select(-Direction, -Score)

  comparison_result <- full_join(DE_res_scores,
                                 misty_res_scores,
                                 by = c("Metabolite", "Cluster")) %>%
    mutate(Source = case_when(
      paste(Metabolite, Cluster) %in% paste(DE_res_scores$Metabolite, DE_res_scores$Cluster) &
        paste(Metabolite, Cluster) %in% paste(misty_res_scores$Metabolite, misty_res_scores$Cluster) ~ "Both",
      paste(Metabolite, Cluster) %in% paste(DE_res_scores$Metabolite, DE_res_scores$Cluster) ~ "DE only",
      paste(Metabolite, Cluster) %in% paste(misty_res_scores$Metabolite, misty_res_scores$Cluster) ~ "MISTY Only",
      TRUE ~ "Unknown"
    ))

  return(comparison_result)
}

SNR_per_region = function(seurat_obj, region){
  DefaultAssay(seurat_obj) = "Metabo"

  seurat_obj@graphs = list()

  features = rownames(seurat_obj)

  expr_data <- FetchData(seurat_obj, vars = c(features, "region"))
  expr_long <- expr_data %>%
    pivot_longer(cols = all_of(features), names_to = "feature", values_to = "expression") %>%
    dplyr::rename(group = region)

  # Compute mean and variance per group-feature combo
  summary_stats <- expr_long %>%
    dplyr::group_by(feature, group) %>%
    dplyr::summarise(
      mean_expr = mean(expression, na.rm = TRUE),
      var_expr = var(expression, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(SNR = mean_expr / var_expr)

  return(summary_stats)
}


Run_fisher_SNR_importance = function(reg_data){
  region_name = unique(reg_data$group)
  fisher_data = reg_data %>%
    dplyr::mutate(High_Importance = ifelse(reg_data$Score >= 2, "Yes", "No"),
                  High_SNR = ifelse(reg_data$SNR >= 2, "Yes", "No")) %>%
    dplyr::select(High_Importance, High_SNR) %>%
    table() %>%
    as.data.frame() %>%
    dplyr::mutate(High_Importance = paste0("High.Imp_", High_Importance),
                  High_SNR = paste0("High.SNR_", High_SNR)) %>%
    pivot_wider(names_from = High_Importance, values_from = Freq) %>%
    column_to_rownames("High_SNR")
  fisher_data = fisher_data[c("High.SNR_Yes", "High.SNR_No"),
                            c("High.Imp_Yes", "High.Imp_No")]

  expected = ((fisher_data[1,1] + fisher_data[1,2]) * (fisher_data[1,1] + fisher_data[2,1])) /
    sum(fisher_data)

  FE = fisher_data[1,1] / expected

  fisher_res = fisher.test(fisher_data, alternative = "greater")
  return(data.frame(region = region_name,
                    fisher_pval_log = -log10(fisher_res$p.value),
                    FE_score = FE))
}
