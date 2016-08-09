setwd("~/Projects/DO")
#.libPaths("/home/xchen/DO/scripts/library")
load("results/8-05-16/Svenson_liver_eQTL_scan1.RData", verbose=T)
load("data/DO478_ExprData_4eQTLMapping.Rdata")
load("data/Svenson_4eQTLMapping.Rdata")
source("scripts/my_funtions.R")
library(qtl2scan)
library(qtl2convert)
library(qtl2geno)
library(qtl2plot)
library(intermediate)

#change annot names for mediation scan
names(annotations.rna.478)[4]="pos"
annot=annotations.rna.478[annotations.rna.478$EnsemblID%in%colnames(expr),]

####################find significant correlations between expression and pc1 phenptypes#########################
pc1=read.csv("results/7-20-16/pc1.csv")
pc2=read.csv("results/7-20-16/pc2.csv")

###principle component 1###
sigcor_pc1=data.frame()
for(i in 7:ncol(phenotype[,as.character(pc1$X)])){
  cor=mycorr(expr, phenotype[,as.character(pc1$X)][,i])
  cor=cor[which(abs(cor)>0.5)]
  geneName=annot$Gene[as.character(annot$EnsemblID)%in%as.character(names(cor))]
  name=rep(colnames(phenotype)[i], length(cor))
  geneID=names(cor)
  sigcor=rbind(sigcor, cbind(geneID, geneName, name, cor))
  print(i)
  print(paste("length =" ,length(cor)))
}
sigcor$cor=as.numeric(as.character(sigcor$cor))
sigcor=sigcor[order(abs(sigcor$cor),  decreasing =T),]
#genes that are signoficantly correltated with PC1
sort(table(sigcor$geneName))


###principle component 2##
sigcor_pc2=data.frame()
for(i in 7:ncol(phenotype[,as.character(pc2$X)])){
  cor=mycorr(expr, phenotype[,as.character(pc2$X)][,i])
  cor=cor[which(abs(cor)>0.5)]
  geneName=annot$Gene[as.character(annot$EnsemblID)%in%as.character(names(cor))]
  name=rep(colnames(phenotype)[i], length(cor))
  geneID=names(cor)
  sigcor=rbind(sigcor, cbind(geneID, geneName, name, cor))
  print(i)
  print(paste("length =" ,length(cor)))
}
sigcor$cor=as.numeric(as.character(sigcor$cor))
sigcor=sigcor[order(abs(sigcor$cor),  decreasing =T),]
#genes that are signoficantly correltated with PC2
sort(table(sigcor$geneName)) 

################done with pcas for now###########################

#####load in phenotype scan1 data#########
load("results/7-26-16/Svenson_pheno_add.RData", verbose=T)
load("results/7-26-16/Svenson_pheno_intsex.RData", verbose=T)
load("results/7-26-16/Svenson_pheno_intdiet.RData", verbose=T)
###########
#get the difference between the lod peaks
pheno.scan1.dietdiff=pheno.scan1.add
pheno.scan1.sexdiff=pheno.scan1.add
pheno.scan1.dietdiff$lod=pheno.scan1.intdiet$lod-pheno.scan1.add$lod
pheno.scan1.sexdiff$lod=pheno.scan1.intsex$lod-pheno.scan1.add$lod
##

##find the peaks with a max diff of 6 or more
index_scan1=c()
maxlod=data.frame()
for(i in 1:ncol(pheno.scan1.dietdiff$lod)){
  for(j in 1:20){
    max=max_scan1(pheno.scan1.dietdiff, lodcolumn=i, chr=j)
    max$pheno=names(max)[3]
    max$marker=rownames(max)
    names(max)[3]="LOD_scores"
    if(max$LOD_scores>=6){
      maxlod=rbind(maxlod, max)
      index_scan1=c(index_scan1,i)
      print(max)
      }
  print(i)
  #print(dim(maxlod))
  }
}
maxlod.dietdiff=maxlod


##find the peaks with a max diff of 6 or more
index_scan1=c()
maxlod=data.frame()
for(i in 1:ncol(pheno.scan1.sexdiff$lod)){
  for(j in 1:20){
    max=max_scan1(pheno.scan1.sexdiff, lodcolumn=i, chr=j)
    max$pheno=names(max)[3]
    max$marker=rownames(max)
    names(max)[3]="LOD_scores"
    if(max$LOD_scores>=6){
      maxlod=rbind(maxlod, max)
      index_scan1=c(index_scan1,i)
      print(max)
      print(j)
    }
    print(i)
    #print(dim(maxlod))
  }
}
maxlod.sexdiff=maxlod

#find the common phenotypes that have both significant sex and diet interactions
common_sex=maxlod.sexdiff[maxlod.sexdiff$pheno%in%maxlod.dietdiff$pheno,]
common_diet=maxlod.dietdiff[maxlod.dietdiff$pheno%in%maxlod.sexdiff$pheno,]
common_diet[common_sex$chr==common_diet$chr,]
common_sex[common_sex$chr==common_diet$chr,]

########USE UPDATED FUNCTION FOR MEDIATION SCAN WITH INTERACTIVE COVARIATE###############
#do mediation scans on with diet interactive covariate
meddiet=vector("list")
for(i in 1:nrow(common_diet)){
  meddiet[[common_diet$pheno[i]]]=mediation.scan(target =  phenotype[,common_diet$pheno[i]],
                     mediator = expr,
                     annotation = annot,
                     covar = covar,
                     int="diet",
                     qtl.geno = probs_doqtl[,,common_diet$marker[i]])
}
#do mediation scans on with sex interactive covariate
medsex=vector("list")
for(i in 1:nrow(common_sex)){
  medsex[[common_sex$pheno[i]]]=mediation.scan(target =  phenotype[,common_sex$pheno[i]],
                                                     mediator = expr,
                                                     annotation = annot,
                                                     covar = covar,
                                                     int="sex",
                                                     qtl.geno = probs_doqtl[,,common_sex$marker[i]])
}

#plot out as pdf
pdf("results/scan1diff_medscan.pdf")
par(mfrow=c(2,2))
for(i in 1:nrow(common_diet)){
  plot(pheno.scan1.sexdiff, lodcolumn = names(medsex)[i], 
       main=paste(names(medsex)[i],"intsex - add"), xlab="")
  plot(pheno.scan1.dietdiff, lodcolumn = names(medsex)[i], 
       main=paste(names(medsex)[i],"intdiet - add"), xlab="", ylab="")
  plot(medsex[[i]], main=paste("sexint on", common_sex$marker[i]))
  plot(meddiet[[i]], main=paste("dietint on",common_diet$marker[i]))
}
graphics.off()
  