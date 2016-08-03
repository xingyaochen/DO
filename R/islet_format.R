#script to convert the islet RNA seq data into qtl2 format
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
print(args)
directory = args[1]

#a directory that contains the derived RNASeq data
seqPath = args[2]
#name of files that contains genotyping probabilities
genoprobName = args[3]
#name of files that contains RNASeq total counts
exprName = args[4]
#path to phenotyping data
phenoPath = args[5]
#name of csv that contains covariate data
covarName = args[6]
#path to export the final data (in .RData format)
outPath = args[7]

#get rz.transform function
rz.transform <- function(y) {
  rankY = rank(y, ties.method = "average", na.last = "keep")
  rzT = qnorm(rankY / (length(na.exclude(rankY)) + 1))
  rzT
}
#make function to get snp64K marker data
getSNP = function() {
  if (!require("RCurl"))
    install.packages("RCurl")
  if (!require("foreign"))
    install.packages("foreign")
  library(RCurl)
  library(foreign)
  url = "https://raw.githubusercontent.com/16xchen/DO/master/data/snps64K.csv"
  data = getURL(url, .opts = list(ssl.verifypeer = FALSE))
  return(read.csv(textConnection(data)))
}



.libPaths(paste(directory, "scripts/library", sep = "/"))
if (!require("abind")) {
  install.packages("abind",  repos = 'http://cran.rstudio.com/')
}

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



#retrieve snp marker info from my github
snps = getSNP()
#convert chrX to chr20 for calculation purposes
snps$Chr = gsub("X", "20", snps$Chr)

fileList <- list.files(path = seqPath)
if (!is.null(fileList)) {
  tab1 <-
    read.table(
      paste(seqPath, fileList[1], genoprobName, sep = "/"),
      sep = "\t",
      header = FALSE
    )
  tab2 <-
    read.table(
      paste(seqPath, fileList[2], genoprobName, sep = "/"),
      sep = "\t",
      header = FALSE
    )
  probs = abind(tab1, tab2, along = 3)
  
  for (i in fileList[3:length(fileList)]) {
    #load in probs data
    probs.mouse <-
      read.table(paste(seqPath, i, genoprobName, sep = "/"),
                 sep = "\t",
                 header = FALSE)
    #set colnames
    colnames(probs.mouse) <-
      c("A", "B", "C", "D", "E", "F", "G", "H")
    #bind the data to the already existing data
    probs <- abind(probs, probs.mouse, along = 3)
    print(i)
    print(dim(probs))
  }
  
  expr <- data.frame() #just get column 10
  
  for (i in fileList) {
    #load in expression data
    expr.mouse <-
      read.table(paste(seqPath, i, exprName, sep = "/"),
                 sep = "\t",
                 header = T)
    #set expression data rownames to ensembl ids
    #transpose the data so it's a single row and markers are in the columns
    totalexpr.mouse <- as.numeric(as.character(expr.mouse[, 10]))
    #bind the data together
    expr = rbind(expr, t(totalexpr.mouse))
    print(i)
    print(dim(expr))
  }
  colnames(expr) = expr.mouse$locus
  #set mouse id to be y axis, strains to be x, and markers to be z
  probs <- aperm(probs, c(3, 2, 1))
  #identify mice in both sets by the folder name
  rownames(probs) <- fileList
  rownames(expr) <- fileList
  dimnames(probs)[[3]] = snps$SNP_ID
}
#rz transform expression data
expr = as.matrix(expr)
for (i in 1:ncol(expr)) {
  expr[, i] = rz.transform(expr[, i])
}
probs_doqtl=probs
#set working directory
setwd(directory)

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
annot = annot[match(names(expr), annot$ensembl_gene_id), ]

#this is the information on each sample (sex, diet, gen)
cov = read.csv(covarName, header = T)
#match order of samples to expression data
cov = cov[match(rownames(probs), covar[, 1])]

covar = data.frame(ncol = 3)
for (i in 1:ncol(cov)) {
  c = as.matrix(as.numeric(cov[, i]) - 1)
  covar = cbind(covar, c)
}
covar = covar[, -1]
covar = as.matrix(covar)
colnames(covar) = colnames(cov)
rownames(covar) = rownames(cov)


#load in phenotype rdata
ob = load(phenoPath, verbo = T)
#rename phenotype data
phenotype = ob[1]


#convert probs to qtl2 formatted probs
probs = probs_doqtl_to_qtl2(
  probs_doqtl,
  map = snps,
  chr_column = "Chr",
  pos_column = "cM",
  marker_column = "SNP_ID"
)

#get kinship and sex covars
k <- calc_kinship(probs_qtl2, "loco", cores = 16)
sex <- (covar$Sex == "M") * 1
names(sex) = covar[, 1]

#save all the datasets into an RData file
save(snps, annot, covar, expr, probs, k, probs, probs_doqtl, phenotype, file = outPath)
