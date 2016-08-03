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

#set working directroy
setwd(directory)
###########begin making dataset for eQTL viewer########################

###Coef Scans#########
if (!is.null(probs) && !is.null(snps) && !is.null(expr)) {
  #do the coef thing
  scanc = vector("list")
  for (j in 1:ncol(expr)) {
    allmarker = data.frame()
    for (i in 1:20) {
      out2 = scan1coef(probs[, as.character(i)],
                       expr[, j],
                       k[[i]],
                       covar)
      allmarker <- rbind(allmarker, out2$coef)
      print(paste("coef scans on chr", i, " of expression", j))
    }
    scanc[[j]] = allmarker
  }
  save(scanc, file = outPath)
  #format coefficients into 3D array
  coef.array = abind(scanc[[1]], scanc[[2]], along = 3)
  for (i in 3:length(scanc)) {
    coef.array = abind(coef.array, scanc[[i]])
    print(paste(i, "scans abinded"))
  }
  save(scanc, coef.array, file = outPath)
  #aperm to 3D transpose the coef dataset (has not yet worked)
  coef.array.t = aperm(coef.array, c(3, 1, 2))
  rownames(coef.array.t) = colnames(expr)
  save(coef.array.t, scanc, file = outPath)
}
