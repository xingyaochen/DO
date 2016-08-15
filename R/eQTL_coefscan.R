#scripts to do coefficient scans on expression data
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
print(args)
#directory
directory = args[1]
#Path to the qtl2 formatted data
filePath = args[2]
#name of the qtl2 coefscan results ad an RData object
outFile = args[3]
from=args[4]
to=args[5]
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

load(filePath)
#print(expr)
#set working directroy
setwd(directory)
###########begin making dataset for eQTL viewer########################
print(class(expr))
###Coef Scans#########
if (!is.null(probs) && !is.null(snps) && !is.null(expr)) {
  #do the coef thing
for(j in from:to){
    coef_scan = data.frame()
    for (i in 1:20) {
      out2 = scan1coef(probs[, as.character(i)],
                       expr[, j],
                       k[[i]],
                       covar)
      		coef_scan <- rbind(coef_scan, out2$coef)
      		print(paste("coef scans on chr", i, " of expression",j))
		    }
	print(j)
	coef_scan=t(coef_scan)   
	 save(coef_scan, file = paste(outFile, j, ".RData", sep=""))
	 }
}
