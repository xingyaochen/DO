.libPaths("/home/xchen/DO/scripts/library")
load("../../hpcdata/gac/derived/Svenson_DO850/DO478_ExprData_4eQTLMapping.Rdata")
load("DO/data/Svenson_4eQTLMapping.Rdata")

source("https://bioc.ism.ac.jp/biocLite.R")
biocLite("DOQTL")
library(DOQTL)

##########
#change probs name and calculate kinship for DOQTL
probs=probs.478
K = kinship.probs(probs, snps = snps, bychr = TRUE)

#make interactive covars
sexcovar=as.matrix(covar[,1])
colnames(sexcovar)=colnames(covar)[1]
dietcovar=as.matrix(covar[,4])
colnames(dietcovar)=colnames(covar)[4]

#do scan one 
#addtive
pheno.scanone.add = scanone(pheno = phenotype, pheno.col = 7:ncol(phenotype), probs = probs, K = K, 
             addcovar = covar, snps = snps)

#intsex
pheno.scanone.intsex = scanone(pheno = phenotype, pheno.col = 7:ncol(phenotype), probs = probs, K = K,
              addcovar = covar,intcovar=sexcovar, snps = snps)
save(pheno.scanone.add,pheno.scanone.intsex,file="DO/results/Svenson_doqlt_scans_2.RData")

#intdiet
pheno.scanone.intdiet = scanone(pheno = phenotype, pheno.col = 7:ncol(phenotype), probs = probs, K = K,
              addcovar = covar[,-4], intcovar=dietcovar, snps = snps)

save(pheno.scanone.add,pheno.scanone.intdiet, pheno.scanone.intsex, file="DO/results/Svenson_doqlt_scans_2.RData")



