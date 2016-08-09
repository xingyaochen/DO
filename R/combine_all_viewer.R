#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
print(args)
#path to your working directory
directory = args[1]
#path to data that was uses for eQTL mapping
dataPath = args[2]
#path to scan1 data
scan1Path = args[3]
#path to coef scan data
coefPath = args[4]
#the name of final data set the gets exported out
outDir = args[5]
#name of the expriment
nameOfExperiment = args[6]

#dataPath="data/Svenson_4eQTLMapping.Rdata"
#scan1Path="results/liver_DO_Svenson/Svenson_eQTL2.RData"
#coefPath="results/Svenson_eQTL2_coeff.RData"
#outPath="results/heart_eQTL_viewer.RData"

.libPaths(paste(directory, "scripts/library", sep = "/"))
if (!require("abind")) {
  install.packages("abind",  repos = 'http://cran.rstudio.com/')
}

if (!require(c("devtools", "RcppEigen"))) {
  install.packages(c("devtools", "RcppEigen"),
                   repos = 'http://cran.rstudio.com/')
}

library(devtools)

if (!require("biomaRt")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("biomaRt")
}


library(devtools)

install_github(paste0("rqtl/qtl2", c("geno",
                                     "scan", "plot", "convert")))

library(abind)
library(biomaRt)
library(qtl2geno)
library(qtl2scan)
library(qtl2plot)
library(qtl2convert)



#scripts for bringing all things together into eQTL viewer dataset
setwd(directory)
load(dataPath)

load(scan1Path, verbose = T)
#print(class(scan1))
#load(coefPath, verbose = T)
load(dataPath, verbose = T)


#retrieve gene annotations from ensembl
ensembl = useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
#get the dataset
annot <- getBM(
  attributes = c(
    'ensembl_gene_id',
    'chromosome_name',
    'start_position',
    'end_position',
    'gene_biotype',
    'external_gene_name',
    'description'
  ),
  filters = 'chromosome_name',
  values = c(1:19, 'X'),
  mart = ensembl
)

#get only the ones with matching expression data
#annot = annot[match(colnames(expr), annot$ensembl_gene_id),]
annot = annot[annot$ensembl_gene_id%in%colnames(expr),]


dataset = list()

###features###########
#annotations are the gene annotatons with ensemblID, chr, location, gene description,
#and gene name
features = data.frame(
  feature_id = annot$ensembl_gene_id,
  group_id = annot$ensembl_gene_id,
  chrom = annot$chromosome_name,
  location = annot$start_position,
  name = annot$external_gene_name,
  description = annot$description
)
dataset$features = features
#print(head(dataset$features))



###markers###########
#snps are set as the  snp info
markers = data.frame(
  marker_id = snps$SNP_ID,
  chrom = snps$Chr,
  location = snps$Mb_NCBI38
)
#put everything into dataset
dataset$markers = markers
#print(head(dataset$markers))


lod = list()
###QTL scans##########
lod$lod = liver_eQTL_scan1$lod
dataset$lod = lod
#print(head(dataset$lod))
coef = list()
strains = data.frame(strain_id = probs$alleles, name = names(CCcolors))
coef$strains = strains
#coef$coef = coef.array.t
#insert coef into dataset
dataset$coef = coef
print(attributes(dataset$coef))

###samples############
samples = data.frame(sample_id = rownames(phenotype),
                     name = rownames(phenotype))
dataset$samples = samples
#print(dataset$samples)

###phenotypes#########
if (!is.null(phenotype_descr)) {
  factors <- data.frame(
    factor_id = phenotype_descr$name,
    name = phenotype_descr$name,
    description = phenotype_descr$description
  )
  dataset$phenotypes$factors = factors
}
if (!is.null(phenotype)) {
  dataset$phenotypes$phenotypes = phenotype
}
#print(head(dataset$phenotype))

###genotypes##########

###expression#########
dataset$expression$expression <- t(expr)
#print(head(expr))

#make attributes
attributes <- list()
attributes = attributes(dataset)
attributes$nameOfExperiment = nameOfExperiment
attributes$DfA <- 7
attributes$DfX <- 14
#
attributes(dataset) <- attributes
print(attributes(dataset))
###save data#########


setwd(outDir)
head(dataset$features)
write.table(dataset$features, file="features.txt")


print(head(dataset$markers))
write.table(dataset$markers, file="markers.txt")
system("cp markers.txt markers-compressed.txt")
system("gzip markers-compressed.txt")
system("rm markers.txt")


dim(dataset$lod$lod)
write.table(dataset$lod$lod, file="lod/lod.txt")
system("gzip lod/lod.txt")


print(head(dataset$expression$expression))
write.table(dataset$expression$expression, file="expression/expression.txt")
system("gzip expression/expression.txt")


print(head(dataset$coef$strains))
write.table(dataset$coef$strains, file="coef/strains.txt")

print(head(dataset$phenotype$factors))
write.table(dataset$phenotypes$factors, file="phenotypes/factors.txt")

head(dataset$phenotypes$phenotypes)
write.table(dataset$phenotypes$phenotypes, file="phenotypes/phenotypes.txt")





