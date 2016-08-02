.libPaths("/home/xchen/DO/scripts/library")

#install packages if not already installed

if(!require(devtools))
  install.packages(c("devtools", "yaml",
                     "jsonlite", "data.table", "RcppEigen"), 
                   repos='http://cran.rstudio.com/')
library(devtools)


if(!require("qtl2geno")|!require("qtl2scan")|!require("qtl2plot")|
   !require("qtl2convert"))
  install_github(paste0("rqtl/qtl2", c("geno",
                                       "scan", "plot", "convert")))

#library(biomaRt)
library(qtl2scan)
library(qtl2plot)
library(qtl2convert)



setwd("/home/xchen/DO")

load("data/Svenson_4eQTLMapping.Rdata")
load("results/liver_DO_Svenson/Svenson_pheno_add.RData")
load("results/liver_DO_Svenson/Svenson_pheno_intsex.RData")
load("results/liver_DO_Svenson/Svenson_pheno_intdiet.RData")




sexcovar=as.matrix(covar[,1])
colnames(sexcovar)=colnames(covar)[1]

dietcovar=as.matrix(covar[,4])
colnames(dietcovar)=colnames(covar)[4]



#initialize empty list of coefficient scans
coefscan_add=vector("list")
coef.scan.each=vector("list")
i=0
if(i>1){
#open pdf stream
pdf("results/Svenson_allchr_coef_add.pdf")
#forloop through every phenotype
for(j in colnames(pheno.scan1.add$lod)){
  #forloop through all 20 chromosomes
  for(i in c(1:20)){
    chr=paste("chr",i,sep="")
    pp=as.matrix(phenotype[,j])
    rownames(pp)=rownames(phenotype)
    colnames(pp)=names(phenotype)[j]
    #do scan1coef
    out2=scan1coef(probs[,as.character(i)],
                   pp,
                   k[[i]],
                   covar)
    plot_coefCC(out2, main=paste(j, "effects on Chr",i, "| Addcovar"), scan1_output = pheno.scan1.add[,j],xlab="")
    #rbinding all the chromosomes together 
    print(paste( "coef scans on chr",i," of expression",j))
    coef.scan.each[[chr]]=out2
  }
  coefscan_add[[j]]=coef.scan.each
  #then put all the rbinded scans into the large list
}
#turn off pdf stream
graphics.off()
save(coefscan_add, file="Svenson_allcoef.RData")
}



#initialize empty list of coefficient scans
coefscan_intsex=vector("list")
coef.scan.each=vector("list")

#open pdf stream
pdf("results/Svenson_allchr_coef_intsex.pdf")
#forloop through every phenotype
for(j in colnames(pheno.scan1.add$lod)){
  #forloop through all 20 chromosomes
  for(i in c(1:20)){
    chr=paste("chr",i,sep="")
    pp=as.matrix(phenotype[,j])
    rownames(pp)=rownames(phenotype)
    colnames(pp)=names(phenotype)[j]
    #do scan1coef
    out2=scan1coef(probs[,as.character(i)],
                   pp,
                   k[[i]],
                  addcovar=covar[,-1],
		intcovar=sexcovar)
    plot_coefCC(out2, main=paste(j, "effects on Chr",i,"| Intcovar=Sex"), scan1_output = pheno.scan1.intsex[,j],xlab="")
    #rbinding all the chromosomes together 
    print(paste( "coef scans on chr",i," of expression",j))
    coef.scan.each[[chr]]=out2
  }
  coefscan_intsex[[j]]=coef.scan.each
  #then put all the rbinded scans into the large list
}
#turn off pdf stream
graphics.off()
save(coefscan_add, coefscan_intsex, file="Svenson_allcoef.RData")



#initialize empty list of coefficient scans
coefscan_intdiet=vector("list")
coef.scan.each=vector("list")

#open pdf stream
pdf("results/Svenson_allchr_coef_intdiet.pdf")
#forloop through every phenotype
for(j in colnames(pheno.scan1.add$lod)){
  #forloop through all 20 chromosomes
  for(i in c(1:20)){
    chr=paste("chr",i,sep="")
    pp=as.matrix(phenotype[,j])
    rownames(pp)=rownames(phenotype)
    colnames(pp)=names(phenotype)[j]
    #do scan1coef
    out2=scan1coef(probs[,as.character(i)],
                   pp,
                   k[[i]],
                  addcovar=covar[,-4],
                intcovar=dietcovar)
    plot_coefCC(out2, main=paste(j, "effects on Chr",i,"| Intcovar=Diet"), scan1_output = pheno.scan1.intdiet[,j], xlab="")
    #rbinding all the chromosomes together
    print(paste( "coef scans on chr",i," of expression",j))
    coef.scan.each[[chr]]=out2
  }
  coefscan_intsex[[j]]=coef.scan.each
  #then put all the rbinded scans into the large list
}
#turn off pdf stream
graphics.off()
save(coefscan_add, coefscan_intsex, coefscan_intdiet,file="results/Svenson_allcoef.RData")



