setwd("/home/xchen/DO")
.libPaths("/home/xchen/DO/scripts/library")
load("results/liver_DO_Svenson/Svenson_liver_eQTL_scan1.RData")
load("data/Svenson_4eQTLMapping.Rdata")
library(qtl2scan)
index_scan1=c()
maxlod=data.frame()
for(i in 1:ncol(liver_eQTL_scan1$lod)){
  max=max_scan1(liver_eQTL_scan1, lodcolumn=i, chr=9)
  name=annot$Gene[which(as.character(annot$EnsemblID)==as.character(names(max)[3]))]
  max$gene_id=names(max)[3]
  max$gene_name=name
  max$marker=rownames(max)
  names(max)[3]="LOD_scores"
  if(max$LOD_scores>=10){
print(max)
    maxlod=rbind(maxlod, max)
    index_scan1=c(index_scan1,i)
  }
  print(i)
  print(dim(maxlod))
}
dim(maxlod)
maxlod_eQTL=maxlod
save(maxlod_eQTL, file="results/Svenson_eQTL_maxlod.RData")
