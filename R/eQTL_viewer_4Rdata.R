#Xingyao Chen + Clifton Jeffery
#7/05/16
#load("../../hpcdata/gac/derived/JAC_DO_Heart_RNASeq/emase_m4_gbrs/rdata/DO192_heart_emase_m4.RData")

#data formatting for inputting into eQTL viewer
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args)
#path to your working directory
directory=args[1]
#names of the objects in your data
#name of phenotype data frame
pheno=args[2]
#names of mRNA expression dataframe
expr=args[3]
#name of genotype probability dataframe
probs=args[4]
#name of SNP dataframe
snps=args[5]
#path to the data directory
filePath=args[6]
nameOfExperiment=args[7]




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


#snps$Chr=gsub("X","20", snps$Chr)

#retrieve snp marker info from my github



#specify the location to which libraries are installed (make a folder named library in 
#your working directory)

#install packages if not already installed

.libPaths("/home/xchen/DO/scripts/library")
#install.packages(c("devtools", "yaml",
 #                 "jsonlite", "data.table", "RcppEigen"), repos='http://cran.rstudio.com/')
install.packages("abind", repos='http://cran.rstudio.com/')
library(abind)
library(devtools)
install_github(paste0("rqtl/qtl2", c("geno", 
                                     "scan", "plot", "convert")))


if(!require("biomaRt")){
  source("https://bioconductor.org/biocLite.R")
  biocLite("biomaRt")
}



library(qtl2geno)
library(qtl2scan)
library(qtl2plot)
library(qtl2convert)
library(abind)
library(biomaRt)


#load in the RData for filepath
load("../../hpcdata/gac/derived/JAC_DO_Heart_RNASeq/emase_m4_gbrs/rdata/DO192_heart_emase_m4.RData")
#set current working directory
setwd(directory)

snps=getSNP()
#check to see if sample numbers match
stopifnot(rownames(expr.mrna)==rownames(pheno)|
            rownames(expr.mrna)==rownames(probs))

snps$Chr=gsub("X","20", snps$Chr)

expr=expr.mrna
ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
#get the dataset for gene annotations
annot <- getBM(attributes=c('ensembl_gene_id',
                            'chromosome_name',
                            'start_position','end_position', 'gene_biotype',
                            'external_gene_name','description'),
               filters='chromosome_name', values=c(1:19,'X'),
               mart = ensembl)

#tolowercase phenortype names
names(pheno)=tolower(names(pheno))

#match annotations to gene expression data
annot=annot[match(colnames(expr), annot$ensembl_gene_id),]


###features###########
if(!is.null(annot)){
  dataset=list()
  #annotations are the gene annotatons with ensemblID, chr, location, gene description,
  #and gene name
  features=data.frame(feature_id=annot$ensembl_gene_id,
                      group_id=annot$ensembl_gene_id,
                      chrom=annot$chromosome_name,
                      location=annot$start_position,
                      name=annot$external_gene_name,
                      description=annot$description)
  dataset$features=features
}


###markers###########
if(!is.null(snps)){
  #snps are set as the  snp info
  markers=data.frame(marker_id=snps$SNP_ID, chrom=snps$Chr,
                     location=snps$Mb_NCBI38)
  #put everything into dataset
  
  dataset$markers=markers
}

###QTL scans##########
if(!is.null(probs) && !is.null(snps) && !is.null(expr)){
  if (class(probs) != "calc_genoprob"){
    probs=probs_doqtl_to_qtl2(probs, map=snps,
                              chr_column = "Chr", pos_column = "cM",
                              marker_column = "SNP_ID")
  }
  qtl2apr = genoprob_to_alleleprob(probs)
  
  #get kinship and sex covars
  k <- calc_kinship(probs, "loco", cores=4)
  sex <- (pheno$sex=="M")*1
  names(sex)=rownames(pheno)
print(attributes(expr))
 # print("scan1-ing")  
#   out <- scan1(qtl2apr, expr[,1], k, sex)
lod=list() 
# lod$lod=out$lod
 #   print("finished scan1 expr")
     # save(out,file="results/scan1_eQTL_heart_x-20.RData")
 # dataset$lod=lod
}


###Coef Scans#########
if(!is.null(probs) && !is.null(snps) && !is.null(expr)){
  #do the coef thing
  scanc=vector("list")
  for(j in 1:ncol(expr)){
    allmarker=data.frame()
    for(i in 1:20){
      out2=scan1coef(qtl2apr[,as.character(i)],
                     expr[,j],
                     k[[i]],
                     sex)
      allmarker <- rbind(allmarker, out2$coef)
      print(paste( "coef scans on chr",i," of expression",j))
    }
    scanc[[j]]=allmarker
  }
save(scanc, file="../results/coeff_eQTL_heart_x-20_list.RData")
  #format coefficients into 3D array
  coef.array=abind(scanc[[1]], scanc[[2]], along=3)
  for(i in 3:length(scanc)){
print(i)   
 coef.array=abind(coef.array, scanc[[i]])
  }
  save(coef.array, file="../results/coeff_eQTL_heart_x-20.RData")

print("aperm-ing")
  coef.array.t=aperm(coef.array, c(3,1,2))
  rownames(coef.array.t)=colnames(expr)
print("done aperms")
  save(coef.array.t, file="../results/coeff_eQTL_heart_x-20_aperms.RData")
  
  #make the coef list
  coef=list()
  strains=data.frame(strain_id=probs$alleles, name=names(CCcolors))
  coef$strains=strains
  coef$coef=coef.array.t
  #insert coef into dataset
  dataset$coef=coef
}
if(0>1){
###Attributes#########
if(!is.null(nameOfExperiment)){
  attributes <- data.frame()
  dataset <- data.frame()
  
  DfA <- 7
  DfX <- 14
  
  attributes <- rbind(attributes, nameOfExperiment)
  attributes <- cbind(attributes, DfA, DfX)
  attributes <- as.list(attributes)
  names(attributes) <- c("name", "DfA", "DfX")
  attributes(dataset) <- attributes
}


###samples############
if(!is.null(phenotype)){
  samples=data.frame(sample_id=rownames(phenotype),
                     name=rownames(phenotype), 
                     description=phenotype$Coat.Color)
  dataset$samples = samples
}


###phenotypes#########
if(!is.null(phenotypeDescriptions) && !is.null(phenotype)){
  factors <- data.frame(factor_id=phenotypeDescriptions$name, 
                        name=phenotypeDescriptions$name,
                        description=phenotypeDescriptions$description)
  dataset$phenotypes$factors = factors
  dataset$phenotypes$phenotypes = newpheno
}

###genotypes##########
if(!is.null(genotype)){
  dataset$genotype$genotype <- t(genotype)
}


###expression#########
if(!is.null(expr.mrna)){
  dataset$expression$expression <- t(expr.mrna)
}

###save data##########
save(dataset, file="results/heart_DO/heart_viewer.Rdata")


}














