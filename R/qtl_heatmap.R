library(DOQTL)
setwd("~/Projects/DO")
#load in the phenotype scanone data
load("results/7-26-16/Svenson_pheno_add.RData")
load("results/7-26-16/Svenson_pheno_intdiet.RData")
load("results/7-26-16/Svenson_pheno_intsex.RData")
load("data/Svenson_4eQTLMapping.Rdata")

#read the phenotypes within Principle Component 1 and 2 
pc1=read.csv("results/7-20-16/pc1.csv")
pc2=read.csv("results/7-20-16/pc2.csv")

#make the titles for the plots
PC=c(rep("PC1",3),rep("PC2",3))
factor=rep(c("Additive","Sex Interactive", "Diet Interactive"),2)
titles=c()
for(i in 1:6){
  titles[i]=paste(PC[i], "QTL Heatmap", "|", factor[i])
}

#Open pdf stream
pdf("results/8-01-16/qtlheatmaps_individualchr_pc1.pdf")
#plot qtl  heatmap of PC1 for each individual chromosome, add, intsex,and intdoet
for(i in 1:20){
  #additive
  ind=colnames(pheno.scan1.add$lod)%in%pc1$X
  lod=cbind(snps[,c(1:2,4)],
            as.data.frame(pheno.scan1.add$lod[,ind]))
  #quartz()
  qtl.heatmap(lod, chr=i)
  if(!is.null(i)){
    title(paste("Chr", i, titles[1]))
  }else{
    title(titles[1])
  }
  
  #intsex
  ind=colnames(pheno.scan1.intsex$lod)%in%pc1$X
  
  lod=cbind(snps[,c(1:2,4)],
            as.data.frame(pheno.scan1.intsex$lod[,ind]))
  qtl.heatmap(lod, chr=i)
  
  if(!is.null(chr)){
    title(paste("Chr", i, titles[2]))
  }else
    title(titles[2])
  
  #intdiet
  ind=colnames(pheno.scan1.intdiet$lod)%in%pc1$X
  
  lod=cbind(snps[,c(1:2,4)],
            as.data.frame(pheno.scan1.intdiet$lod[,ind]))
  
  #quartz()
  qtl.heatmap(lod, chr=i)
  if(!is.null(Chr)){
    title(paste("Chr", i, titles[3]))
  }else{
    title(titles[3])
  }
  print(i)
  #plot.new()
}
graphics.off()


#open pdf stream
pdf("results/8-01-16/qtlheatmaps_individualchr_pc2.pdf")
#plot QTL heatmap of PC2 for each individual chromosome, additive, intsex, and intdiet
for(i in 1:20){
  #additive
  ind=colnames(pheno.scan1.add$lod)%in%pc2$X
  lod=cbind(snps[,c(1:2,4)],
            as.data.frame(pheno.scan1.add$lod[,ind]))
  #quartz()
  qtl.heatmap(lod, chr=i)
  if(!is.null(i)){
    title(paste("Chr", i, titles[4]))
  }else{
    title(titles[4])
  }
  
  #intsex
  ind=colnames(pheno.scan1.intsex$lod)%in%pc2$X
  
  lod=cbind(snps[,c(1:2,4)],
            as.data.frame(pheno.scan1.intsex$lod[,ind]))
  qtl.heatmap(lod, chr=i)
  
  if(!is.null(chr)){
    title(paste("Chr", i, titles[5]))
  }else
    title(titles[5])
  
  
  #intdiet
  ind=colnames(pheno.scan1.intdiet$lod)%in%pc2$X
  
  lod=cbind(snps[,c(1:2,4)],
            as.data.frame(pheno.scan1.intdiet$lod[,ind]))
  
  #quartz()
  qtl.heatmap(lod, chr=i)
  if(!is.null(Chr)){
    title(paste("Chr", i, titles[6]))
  }else{
    title(titles[6])
  }
  print(i)
  # plot.new()
}
graphics.off()










