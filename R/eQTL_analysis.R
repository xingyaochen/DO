setwd("~/Projects/DO")
#.libPaths("/home/xchen/DO/scripts/library")
load("results/8-05-16/Svenson_liver_eQTL_scan1.RData", verbose = T)
load("data/DO478_ExprData_4eQTLMapping.Rdata", verbose = T)
load("data/Svenson_4eQTLMapping.Rdata")
source("scripts/my_funtions.R")
library(qtl2scan)
library(qtl2convert)
library(qtl2geno)
library(qtl2plot)
library(intermediate)
source("scripts/my_funtions.R")
#change annot names for mediation scan
names(annotations.rna.478)[4] = "pos"
#annot=annotations.rna.478[annotations.rna.478$EnsemblID%in%colnames(expr),]

####################find significant correlations between expression and pc1 phenptypes#########################
pc1 = read.csv("results/7-20-16/pc1.csv")
pc2 = read.csv("results/7-20-16/pc2.csv")

###principle component 1###
sigcor_pc1 = data.frame()
for (i in 7:ncol(phenotype[, as.character(pc1$X)])) {
  cor = mycorr(expr, phenotype[, as.character(pc1$X)][, i])
  cor = cor[which(abs(cor) > 0.4)]
  geneName = annot$Gene[as.character(annot$EnsemblID) %in% as.character(names(cor))]
  name = rep(colnames(phenotype)[i], length(cor))
  geneID = names(cor)
  sigcor_pc1 = rbind(sigcor_pc1, cbind(geneID, geneName, name, cor))
  print(i)
  print(paste("length =" , length(cor)))
}
sigcor$cor = as.numeric(as.character(sigcor$cor))
sigcor = sigcor[order(abs(sigcor$cor),  decreasing = T), ]
#genes that are signoficantly correltated with PC1
sort(table(sigcor$geneName))

pc$scores = pc$scores[rownames(pc$scores) %in% rownames(expr), ]
expr = expr[rownames(expr) %in% rownames(pc$scores), ]
cor = mycorr(expr, pc$scores[, 1])
cor = cor[which(abs(cor) > 0.3)]
geneName = annot$external_gene_name[as.character(annot$ensembl_gene_id) %in%
                                      as.character(names(cor))]
#name=rep(colnames(phenotype)[i], length(cor))
geneID = names(cor)
sigcor_pc1 = data.frame(geneID, geneName, cor)
sigcor_pc1[order(abs(sigcor_pc1$cor), decreasing = F), ]

###principle component 2##
sigcor_pc2 = data.frame()
for (i in 7:ncol(phenotype[, as.character(pc2$X)])) {
  cor = mycorr(expr, phenotype[, as.character(pc2$X)][, i])
  cor = cor[which(abs(cor) > 0.5)]
  geneName = annot$Gene[as.character(annot$EnsemblID) %in% as.character(names(cor))]
  name = rep(colnames(phenotype)[i], length(cor))
  geneID = names(cor)
  sigcor = rbind(sigcor, cbind(geneID, geneName, name, cor))
  print(i)
  print(paste("length =" , length(cor)))
}
sigcor$cor = as.numeric(as.character(sigcor$cor))
sigcor = sigcor[order(abs(sigcor$cor),  decreasing = T), ]
#genes that are signoficantly correltated with PC2
sort(table(sigcor$geneName))

####find genes with high lod peaks ok chr 9
index = which(annot$chromosome_name == 9)
chr_scan1 = liver_eQTL_scan1[, index]
class(chr_scan1)
dim(chr_scan1$lod)
index_scan1 = c()
maxlod = data.frame()
for (i in 1:ncol(chr9_scan1$lod)) {
  max = max_scan1(chr9_scan1, lodcolumn = i)
  name = annot$Gene[which(as.character(annot$EnsemblID) == as.character(names(max)[3]))]
  max$gene_id = names(max)[3]
  max$gene_name = name
  max$marker = rownames(max)
  names(max)[3] = "LOD_scores"
  if (max$LOD_scores >= 12) {
    maxlod = rbind(maxlod, max)
    index_scan1 = c(index_scan1, i)
    print(max)
  }
  print(i)
  #print(dim(maxlod))
}
dim(maxlod)
maxlod9_eQTL = maxlod

index = which(annot$chromosome_name == 6)
chr_scan1 = liver_eQTL_scan1[, index]
class(chr_scan1)
dim(chr_scan1$lod)
index_scan1 = c()
maxlod = data.frame()
for (i in 1:ncol(chr_scan1$lod)) {
  max = max_scan1(chr_scan1, lodcolumn = i)
  name = annot$Gene[which(as.character(annot$EnsemblID) == as.character(names(max)[3]))]
  max$gene_id = names(max)[3]
  max$gene_name = name
  max$marker = rownames(max)
  names(max)[3] = "LOD_scores"
  if (max$LOD_scores >= 12) {
    maxlod = rbind(maxlod, max)
    index_scan1 = c(index_scan1, i)
    print(max)
  }
  print(i)
  #print(dim(maxlod))
}
maxlod6_eQTL = maxlod
nm = maxlod9_eQTL$gene_id[maxlod_eQTL$gene_id %in% sigcor_pc1$geneID]
chr9Genes = data.frame(annot[annot$ensembl_gene_id %in% nm, 1:6], LOD_scores =
                         maxlod9_eQTL$LOD_scores[maxlod9_eQTL$gene_id %in% nm],
                       on_chr = maxlod9_eQTL$chr[maxlod9_eQTL$gene_id %in% nm])
head(chr9Genes)
##
nm = maxlod6_eQTL$gene_id[maxlod6_eQTL$gene_id %in% sigcor_pc1$geneID]
chr6Genes = annot[annot$ensembl_gene_id %in% nm, ]
chr6Genes = data.frame(annot[annot$ensembl_gene_id %in% nm, 1:6],
                       LOD_scores = maxlod6_eQTL$LOD_scores[maxlod6_eQTL$gene_id %in%
                                                              nm],
                       on_chr = maxlod6_eQTL$chr[maxlod6_eQTL$gene_id %in% nm])
head(chr6Genes)
load("data/Svenson_4eQTLMapping.Rdata")
plot_scan1(liver_eQTL_scan1, lodcolumn = chr9Genes$ensembl_gene_id[1], main=chr9Genes$external_gene_name[1])
plot_scan1(liver_eQTL_scan1, lodcolumn = chr6Genes$ensembl_gene_id[1], main=chr6Genes$external_gene_name[1])
