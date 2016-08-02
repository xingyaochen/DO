#using qtl2 to do genome scans on Svenson phenotype data
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
directory=args[1]
int=args[2]
print(paste("directory is",args[1]))
print(paste("interactive covariate is",args[2]))
print(paste("files are saved to",args[3]))


if(tolower(int)=="sex"){
  num=1
}else if(tolower(int)=="diet"){
  num=4
}else{
  num=NULL
}          
.libPaths(paste0(directory,"/library"))

if(!require(c("devtools", "RcppEigen"))){
  install.packages(c("devtools", "RcppEigen"), 
                   repos='http://cran.rstudio.com/')
}

library(devtools)
if(!require(paste0("rqtl/qtl2", c("geno",
                                  "scan", "plot", "convert")))){
  install_github(paste0("rqtl/qtl2", c("geno",
                                       "scan", "plot", "convert")))
}

library(biomaRt)
library(qtl2scan)
library(qtl2plot)
library(qtl2convert)
library(qtl2geno)

setwd("/home/xchen/DO")
load("data/Svenson_4eQTLMapping.Rdata")
if(!is.null(num)){
  #make interactive covariates
  intcovar=as.matrix(covar[,num])
  colnames(intcovar)=colnames(covar)[num]
  print("scan1-ing")
  pheno.scan1 <- scan1(probs, phenotype[,7:ncol(phenotype)], k, addcovar=covar,intcovar=intcovar)
  print("done scan1")
}else{
  print("scan1-ing")
  pheno.scan1 <- scan1(probs, phenotype[,7:ncol(phenotype)], k, addcovar=covar)
  print("done scan1")
}
save(pheno.scan1, file=outPath)

#get maximum lod score
#load("results/Svenson_pheno_scan1.RData")
index_scan1=c()
maxlod=data.frame()
for(i in 1:ncol(pheno.scan1$lod)){
  max=max_scan1(pheno.scan1, lodcolumn=i)
  max$pheno=names(max)[3]
  names(max)[3]="LOD_scores"
  if(max$LOD_scores>=7){
    maxlod=rbind(maxlod, max)
    index_scan1=c(index_scan1,i)
  }
}
dim(maxlod)
maxlod
save(maxlod, pheno.scan1, file=outPath)

#plot coef scans if no interactive covariate
if(is.null(num)){
  scan_coef=vector("list")
  #make pdf stream
  pdf("Svenson_highLod_coef.pdf")
  for(i in 1:nrow(maxlod)){
    allmarker=data.frame()
    j=maxlod$pheno[i]
    chr=maxlod$chr[i]
    pp=as.matrix(phenotype[,j])
    rownames(pp)=rownames(phenotype)
    colnames(pp)=j
    out2=scan1coef(probs[,chr],
                   pp,
                   k[[as.numeric(chr)]],
                   covar)
    plot_coefCC(out2, main=paste(j, "effects on Chr",chr), scan1_output = pheno.scan1[,index_scan1[i]])
    legend("topright", col=CCcolors, names(CCcolors), ncol=2, lwd=2, bg="gray95")
    allmarker <- rbind(allmarker, out2$coef)
    print(paste( "coef scans on chr",chr," of ",j))
    scan_coef[[j]]=allmarker
  }
  graphics.off()
  save(scan_coef,maxlod,pheno.scan1, file=outPath)
}
