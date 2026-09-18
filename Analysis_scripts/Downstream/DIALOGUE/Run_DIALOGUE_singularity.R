library(DIALOGUE)

# /mnt/library is the host R library bound in by run_dialogue_singularity.sh
if (!requireNamespace("nnls", lib.loc = "/mnt/library", quietly = TRUE)) {
  install.packages("nnls", lib = "/mnt/library")
}
library(nnls, lib.loc = "/mnt/library")

args <- commandArgs(trailingOnly = TRUE)

input = readRDS(args[1])
outdir = args[2]


# Functions ---------------------------------------------------------------
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
  idx<-unique(DIALOGUE:::get.strsplit(names(R$sig[[1]]),".",1))
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


# Main --------------------------------------------------------------------
param <- DIALOGUE::DLG.get.param(k = 10,
                                 results.dir = outdir,
                                 conf = c("cellQ"), spatial.flag = T,
                                 abn.c = 15,
                                 n.genes = 100)
param$averaging.function = colMeans2


D1_res = DIALOGUE:::DIALOGUE1(rA = input,main = "CT_5_run",
                               param = param)
print("Finished Running DIALOGUE 1\n")

D2_res = DIALOGUE:::DIALOGUE2(rA = input,main = "CT_5_run",
                               results.dir = param$results.dir)

print("Finished Running DIALOGUE 2\n")


saveRDS(D1_res, file.path(outdir, paste0("DIALOGUE1", "_CT_5_run.rds")))
saveRDS(D2_res, file.path(outdir, paste0("DIALOGUE2", "_CT_5_run.rds")))

D2_res = readRDS(file.path(outdir, paste0("DIALOGUE2", "_CT_5_run.rds")))
Res = custom_D3_DIALOGUE(rA_list = input, D2_result = D2_res)

print("Finished Running DIALOGUE 3\n")

saveRDS(Res,file.path(outdir, "updated_CT_5_main_out.rds"))

print("All Done :) ")
