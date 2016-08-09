#script to do scan1s on expression data
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
print(args)
#directory
directory = args[1]
#Path to the qtl2 formatted data
filePath = args[2]
#name of the qtl2 scan1 results ad an RData object
outFile = args[3]
from=args[4]
to=args[5]

.libPaths(paste(directory, "scripts/library", sep = "/"))
#if (!require("abind")) {
#  install.packages("abind",  repos = 'http://cran.rstudio.com/')
#}

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
#library(abind)
library(biomaRt)

#load in the data
load(filePath, verbose = T)
#load("DO/data/DO_4eQTLMapping.RData")

#set working directory
setwd(directory)
###QTL scans##########
if (!is.null(probs) && !is.null(snps) && !is.null(expr)) {
  if (class(probs) != "calc_genoprob") {
    probs = probs_doqtl_to_qtl2(
      probs,
      map = snps,
      chr_column = "Chr",
      pos_column = "cM",
      marker_column = "SNP_ID"
    )
  }
  print("scan1-ing")
  out <- scan1(probs, expr[,from:to], k, covar)
  print("finished scan1")
  save(out, file = outFile)
}
