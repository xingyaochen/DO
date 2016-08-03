#Xingyao Chen + Clifton Jeffery
#7/05/16
#data formatting for heart DO data to do scans
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
print(args)
#path to your working directory
directory = args[1]
#path to the 4eQTL mapping data
dataPath = args[2]

.libPaths(paste(directory, "scripts/library", sep = "/"))
if (!require(c("devtools", "RcppEigen"))) {
 install.packages(c("devtools", "RcppEigen"),
                  repos = 'http://cran.rstudio.com/')
}

library(devtools)
if (!require(paste0("rqtl/qtl2", c("geno",
                                   "scan", "plot", "convert")))) {
  install_github(paste0("rqtl/qtl2", c("geno",
                                       "scan", "plot", "convert")))
}

if (!require("biomaRt")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("biomaRt")
}


library(devtools)
library(qtl2geno)
library(qtl2scan)
library(qtl2plot)
library(qtl2convert)
library(abind)
library(biomaRt)


getSNP=function(){
  if(!require("RCurl"))
    install.packages("RCurl")
  if(!require("foreign"))
    install.packages("foreign")
  library(RCurl)
  library(foreign)
  url= "https://raw.githubusercontent.com/16xchen/DO/master/data/snps64K.csv"
  data = getURL(url, .opts = list(ssl.verifypeer = FALSE))
  return(read.csv(textConnection(data)))
}

#load in data
load(dataPath, verbose = T)

print("the printed names should correcspond to your args")

setwd(directory)
snps=getSNP()
snps$Chr=gsub("X","20", snps$Chr)


#names of the objects in your data
#name of phenotype data frame
phenotype = pheno
#names of mRNA expression dataframe
expr = expr.mrna
#name of genotype probability dataframe
probs_doqtl = probs
#path to save the final 4eQTLMapping data to
outName = args[7]

print(phenotype)
sex=as.numeric(factor(phenotype$sex))-1
gen=as.numeric(factor(phenotype$gen))-1
covar = data.frame(sex=sex, gen=gen)
rownames(covar)=rownames(phenotype)




#retrieve gene annotations from ensembl
ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
#get the dataset
annot <- getBM(attributes=c('ensembl_gene_id',
                            'chromosome_name',
                            'start_position','end_position', 'gene_biotype',
                            'external_gene_name','description'),
               filters='chromosome_name', values=c(1:19,'X'),
               mart = ensembl)

#get only the ones with matching expression data
annot=annot[match(colnames(expr), annot$ensembl_gene_id),]

probs=probs_doqtl_to_qtl2(probs_doqtl, map=snps,
                            chr_column = "Chr", pos_column = "cM",
                            marker_column = "SNP_ID")
k <- calc_kinship(probs, "loco", cores=4)

save(phenotype, snps, expr, probs_doqtl, covar,k, probs, annot, file = outName)
