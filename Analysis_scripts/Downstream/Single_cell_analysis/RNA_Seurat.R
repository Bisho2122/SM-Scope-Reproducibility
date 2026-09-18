library(tidyverse)
library(Seurat)


# Load Data ---------------------------------------------------------------
all_counts = readRDS("../../../Data/Figure 1/Baysor_filtered_count_mat.rds")
# Functions ---------------------------------------------------------------
Run_seurat_analysis = function(count_mat, max_ncount_RNA = 7500, min_cells = 1000, min_features = 8,
                               scale_factor = 10000, vst_mean_thresh = 0, vst_std_variance_thresh = 1,
                               pca_comps = 20, clust_resol = 0.2, get_clust_markers = F, select_feats_vst = T,
                               ions_to_exclude = NULL, scale = F, center = F, use_leiden = F){

  analysis <- CreateSeuratObject(counts = count_mat, min.cells = min_cells, min.features = min_features) #filter genes and cells
  #VlnPlot(analysis, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
  analysis<-subset(analysis, subset = nCount_RNA > 25 & nCount_RNA < max_ncount_RNA) #filter out any cells with less than 20 and more than 1000 counts
  analysis <- NormalizeData(analysis, normalization.method = "LogNormalize", scale.factor = scale_factor)
  analysis <- FindVariableFeatures(analysis, selection.method = "vst")

  z = analysis@assays$RNA@meta.features
  if(select_feats_vst){
    var_feats = z[which(log(z$vst.mean) > vst_mean_thresh & z$vst.variance.standardized >= vst_std_variance_thresh),]
    var_feats = rownames(var_feats)
  }
  else{
    var_feats = rownames(analysis)
  }
  if (!is.null(ions_to_exclude)){
    var_feats = var_feats[which(var_feats %nin% ions_to_exclude)]
  }

  if(pca_comps > length(var_feats)){
    dimensions = c(1:length(var_feats))
  }
  else{
    dimensions = c(1:pca_comps)
  }

  analysis <- ScaleData(analysis, do.scale = scale, do.center = center)
  analysis <- RunPCA(analysis, features = var_feats, approx = F)

  #cluster the cells
  analysis <- FindNeighbors(analysis,
                            dims = dimensions)
  if (use_leiden){
    analysis <- FindClusters(analysis, resolution = clust_resol, algorithm = 4)
  }
  else{
    analysis <- FindClusters(analysis, resolution = clust_resol)
  }
  analysis <- RunUMAP(analysis,
                      dims = dimensions)

  if(get_clust_markers){
    analysis.markers <- FindAllMarkers(analysis,test.use = "wilcox", only.pos = FALSE)
    analysis.markers$diff_pct = analysis.markers$pct.1 - analysis.markers$pct.2
  }
  else{
    analysis.markers = NULL
  }
  return(list("Seurat_obj" = analysis,
              "mean_var_feats" = z,
              "analysis_markers" = analysis.markers))
}
# preprocessing counts and metadata ---------------------------------------------------------
counts = all_counts$count_mat
#counts = as.data.frame(counts)
counts = as.matrix(counts)
counts[is.na(counts)] = 0
rNames = counts[,1]
counts = apply(counts[,colnames(counts)[c(2:ncol(counts))]], 2,as.numeric)
rownames(counts) = rNames

# Get_cell_type_info ------------------------------------------------------
zeisel_brain_types = readxl::read_xlsx("../../../Data/General/Zeisel_brain_cell_types.xlsx")
colnames(zeisel_brain_types) = c("clust_id", "clust_name","Description", "Region",
                                 "development_compartment", "Neurotransmitter", "comment",
                                 "taxonomy_group","probable_location")
loom_cell_types = loomR::connect(filename = "../../../Data/General/l5_all.agg.loom", mode = "r+", skip.validate = T)

brain_cell_type_reference = list()
z = loom_cell_types$col.attrs
for (i in 1:265){
  res = data.frame(clust_name = z$ClusterName[i], marker_genes = z$MarkerGenes[i],
                   marker_robustness = z$MarkerRobustness[i], marker_selectivity = z$MarkerSelectivity[i],
                   marker_specificity = z$MarkerSpecificity[i])
  brain_cell_type_reference[[i]] = res
}
brain_cell_type_reference = dplyr::bind_rows(brain_cell_type_reference)
brain_cell_type_reference = brain_cell_type_reference %>% dplyr::left_join(zeisel_brain_types)

all_zeisel_marker_genes = c()
for (i in 1:265){
  g = trimws(unlist(strsplit(brain_cell_type_reference$marker_genes[i], " ")))
  all_zeisel_marker_genes = c(all_zeisel_marker_genes, g)
}
all_zeisel_marker_genes = unique(all_zeisel_marker_genes)

cell2location_markers = list(
  'Neurons' =  c('Rbfox3', "Snap25", "Ndrg4"), # NeuN
  'VGLUT1 excitatory' = 'Slc17a7',
  'VGLUT2 excitatory' = 'Slc17a6',
  'cortical.layers' = c('Layer I inhibitory' = 'Reln',
                        'Layer II/III-IV excitatory.1' = 'Cux1',
                        'Layer II/III-IV excitatory.2' = 'Cux2',
                        'Layer IV excitatory' = 'Rorb',
                        'Layer Va excitatory' = 'Pcp4',
                        'Layer Vb excitatory' = 'Htr2c',
                        'Layer VIa excitatory' = 'Oprk1',
                        'Layer VIb excitatory' = 'Nr4a2',
                        "Bcl11b", "Foxp2",
                        "Claustrum" = "Synpr"),
  'hippocampus' = c('Neurod6', 'Vipf3', 'Ahcyl2', 'Nr3c2'),
  'amygdala' = c('Kcng1', 'Lypd1', 'Ndst4', 'Kcnk2'),
  "thalamus" = c('Ptpn3', 'Ptpn4'),
  'GABAergic' = 'Gad1',
  'Cholinergic' = 'Ache',
  'Interneuron' = c('Sst', 'Htr3a', 'Vip', "Pvalb","Lamp5", 'Meis2'),
  "Nb" = c("Dcx", 'Sox2', "Pax6"),
  'Astrocytes' = c('Gfap', 'Aqp4','Slc1a3', 'Mfge8', 'Agt' ),
  'Oligodendrocytes' = c('Plp1','Mog'),
  'OPC' = c('Pdgfra', "Bmp4"),
  'Endothelial cells' = c('Apold1', "Flt1", 'Tek'),
  'Microglia' = c('Ccl4','Ccl3','Tmem119', "Ptprc")
)
all_cell2location_markers_1 = unlist(cell2location_markers)

Tsaic_2016_celltypes = readxl::read_xlsx("../../../Data/General/Tsaic_2016_cell_types.xlsx")
all_Tsaic_markers = c()
for (i in 1:nrow(Tsaic_2016_celltypes)){
  all_Tsaic_markers = c(all_Tsaic_markers, unlist(strsplit(Tsaic_2016_celltypes$`Present Markers`[i], ",")))
}
all_Tsaic_markers = unique(all_Tsaic_markers)

cell2location_markers_S7 = list(
  "Excitatory|Claustrum" = c("Nr4a2", "Synpr"),
  "Excitatory|Amygdala" = c("Lypd1", "Kcng1", "C1ql3"),
  "Excitatory|hippocampus CA1" =  c("Rprml","Wipf3","Neurod6"),
  "Excitatory|hippocampus CA3" = c("Nr4a3","Nptx1","Rnf182"),
  "Excitatory|thalamus" = c("Synpo2", "Ptpn3","Slc17a6"),
  "Excitatory|hypothalamus" = c("Prkch","Ramp3","Shox2"),
  "Excitatory|cortex" = c("Sox5","Htr2c","Kcnk2","Dkk3","Rorb","Foxp2","Bcl11b","Cux2",
                          "Thsd7a","Cux1"),
  "Inhibitory|habenula" = c("Nwd2","Lrrc55","Syt9"),
  "Inhibitory|cortex/hippocampus" = c("Pvalb","Vip","Lamp5","Dlx6os1","Kcnip1","Sst","Vwc2"),
  "Inhibitory|striatum" = c("Drd1","Rgs9","Nexn"),
  "Inhibitory|thalamus" = c("Esrrg","Ubash3b","Nova1","Scg2","Fign","Syt2","Tmem130"),
  "Microglia" = c("Ptprc"),
  "Neuroblasts" = c("Sox2","Dcx","Pax6"),
  "Oligodendrocytes" = c("Mog","Cnksr3","Plp1"),
  "OPCs" = c("Bmp4","Pdgfra"),
  "Astrocytes" = c("Aldoc","Slc1a3","Hepacam")
)
zeisel_markers = list()
for (i in 1:nrow(brain_cell_type_reference)){
  x = unlist(strsplit(brain_cell_type_reference$marker_genes[i]," "))
  if(brain_cell_type_reference$Description[i] %in% names(zeisel_markers)){
    zeisel_markers[[brain_cell_type_reference$Description[i]]] = unique(c(zeisel_markers[[brain_cell_type_reference$Description[i]]],x))
  }
  zeisel_markers[[brain_cell_type_reference$Description[i]]] = x
}
all_ref_cell_types = unlist(list(cell2location_markers, cell2location_markers_S7,zeisel_markers), recursive = F)

all_cell_type_markers = list()
for (i in 1:length(all_ref_cell_types)){
  res = data.frame(cell_type = names(all_ref_cell_types)[i],
                   marker_genes = paste(all_ref_cell_types[[i]], collapse = ","),
                   num_marker_genes = length(all_ref_cell_types[[i]]))
  all_cell_type_markers[[i]] = res
}
all_cell_type_markers = dplyr::bind_rows(all_cell_type_markers)
all_cell_type_markers = all_cell_type_markers[!is.na(all_cell_type_markers$cell_type),]
write.csv(all_cell_type_markers, "../../../Data/Figure 1/cell_type_markers.csv", row.names = F)

x = intersect(rownames(counts), unlist(all_ref_cell_types))

Resolve_cell_types = list()
for (i in 1:length(all_ref_cell_types)){
  if(length(intersect(all_ref_cell_types[[i]], x)) != 0){
    res = data.frame(cell_type = names(all_ref_cell_types)[i],
                     resolve_marker_genes = paste(intersect(all_ref_cell_types[[i]], x), collapse = ","),
                     num_marker_genes = length(intersect(all_ref_cell_types[[i]], x)))
    Resolve_cell_types[[i]] = res
  }
  else{
    Resolve_cell_types[[i]] = NULL
  }
}
Resolve_cell_types = dplyr::bind_rows(Resolve_cell_types)

#Manually curate in excel to remove duplicates
write.csv(Resolve_cell_types, "../../../Data/Figure 1/Resolve_cell_types.csv", row.names = F)

Resolve_cell_types = read.csv("../../../Data/Figure 1/Resolve_cell_types.csv")

# Run Seurat --------------------------------------------------------------
input<- counts
# row.names(input)<-input$Gene.ID
# input$Gene.ID<-NULL

analysis <- CreateSeuratObject(counts = input, min.cells = 200, min.features = 8) #filter genes and cells
VlnPlot(analysis, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

seurat_transc = Run_seurat_analysis(count_mat = counts, max_ncount_RNA = 1000, min_cells = 200,
                                    min_features = 4, pca_comps = 10, clust_resol = 0.2,scale = T,
                                    center = T, scale_factor = 10000, select_feats_vst = F, get_clust_markers = T)

cluster_annotation = c("0" = "Excit_CTX",
                       "1" = "Inhib_CTX",
                       "2" = "Oligodendrocytes",
                       "3" = "Excit_VGLUT2",
                       "4" = "Astrocytes",
                       "5" = "Thalamic Neurons",
                       "6" = "Excit_HPF",
                       "7" = "Microglia",
                       "8" = "Inhib_Interneurons",
                       "9" = "Astrocytes",
                       "10" = "Ventricular systems Glial cells")

seurat_transc$Seurat_obj = RenameIdents(seurat_transc$Seurat_obj, cluster_annotation)

saveRDS(list("seurat_res" = seurat_transc,
             "params" = list("PCA_comps" = 10, "Scale" = "True", "Center" = "True", "Clustering_resol" = 0.2,
                             "excluded_genes" = NULL, "max_ncount_ions" = 1000)),
        "../../../Data/Figure 1/Final_seurat_transc_updated.rds")
