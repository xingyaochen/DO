#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args)
directory=args[1]
seqPath=args[2]
genoprobName=args[3]
exprName=args[4]
phenoName=args[5]



.libPaths(paste(directory,"library", sep="/"))
if(!require("abind")){
  install.packages("abind",  repos='http://cran.rstudio.com/')
}

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

if(!require("biomaRt")){
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


#make function to get snp64K marker data
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



#retrieve snp marker info from my github
snps=getSNP()


fileList <- list.files(path=seqPath)
if(!is.null(fileList)){
  tab1 <- read.table(paste(seqPath,fileList[1], genoprobName, sep="/"), sep = "\t", header = FALSE)
  tab2 <- read.table(paste(seqPath,fileList[2], genoprobName, sep="/"), sep = "\t", header = FALSE)
  genoprobs=abind(tab1,tab2, along=3)
  
  for(i in fileList[3:length(fileList)]){
    #load in genoprobs data
    genoprobs.mouse <- read.table(paste(seqPath, i, genoprobName, sep="/"), sep = "\t", header = FALSE)
  #set colnames
    colnames(genoprobs.mouse) <- c("A", "B", "C", "D", "E", "F", "G", "H")
  #bind the data to the already existing data
    genoprobs <- abind(genoprobs, genoprobs.mouse, along = 3)
    print(i)
    print(dim(genoprobs))
  }
  
  expr <- data.frame() #just get column 10
  
  for(i in fileList){
    #load in expression data
    expr.mouse <- read.table(paste(seqPath, i, exprName, sep="/"), sep = "\t", header = T)
    #set expression data rownames to ensembl ids
    #transpose the data so it's a single row and markers are in the columns
    totalexpr.mouse <-as.numeric(as.character(expr.mouse[,10]))
    #bind the data together
    expr=rbind(expr, t(totalexpr.mouse))
    print(i)
    print(dim(expr))
  }
  
  colnames(expr)=expr.mouse$locus
  #set mouse id to be y axis, strains to be x, and markers to be z
  genoprobs <- aperm(genoprobs, c(3,2,1))
  #identify mice in both sets by the folder name
  rownames(genoprobs) <- fileList
  rownames(expr) <- fileList
  dimnames(genoprobs)[[3]]=snps$SNP_ID
}
  


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
annot=annot[match(names(expr), annot$ensembl_gene_id),]

#this is the information on each sample (sex, diet, gen)
covar=read.csv(phenoName, header=T)
#match order of samples to expression data
covar=covar[match(rownames(genoprobs),covar[,1]),]




###########begin making dataset for eQTL viewer########################

#set working directory
setwd(directory)
dataset=list()

###features###########
if(!is.null(annot)){
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
if(!is.null(genoprobs) && !is.null(snps) && !is.null(expr)){
  if (class(genoprobs) != "calc_genoprob"){
    probs=probs_doqtl_to_qtl2(genoprobs, map=snps,
                              chr_column = "Chr", pos_column = "cM",
                              marker_column = "SNP_ID")
  }
  qtl2apr = genoprob_to_alleleprob(probs)
  
  #get kinship and sex covars
  k <- calc_kinship(qtl2apr, "loco", cores=16)
  sex <- (covar$Sex=="M")*1
  names(sex)=covar[,1]
  print("scan1-ing")
  lod=list()
  lodvalue=data.frame(nrow=nrow(snps))
	 pheno=as.matrix(expr[,1])
   	 rownames(pheno)=covar[,1]
   	 colnames(pheno)=colnames(expr)[1]
	allscan1=scan1(qtl2apr, pheno, k, sex, cores=4)

  	for(i in 1:ncol(expr)){
  		  pheno=as.matrix(expr[,i])
  		  rownames(pheno)=covar[,1]
  		  colnames(pheno)=colnames(expr)[i]
  		  out <- scan1(qtl2apr, pheno, k, sex, cores=4)
  		  print(paste("finished scan1 expr", i))
  		  lodvalue=cbind(lodvalue, out$lod)
  		  write.csv(lodvalue,"lod_all.csv")
		if(i>1){
			allscan1=cbind.scan1(allscan1, out)
			save(allscan1, "allscan1.RData")
			}  
		}
  lod$lod=lodvalue[,-1]
  dataset$lod=lod
}


###Coef Scans#########
if(!is.null(probs) && !is.null(snps) && !is.null(expr)){
  #do the coef thing
  scanc=vector("list")
  for(j in 1:ncol(expr)){
    allmarker=data.frame()
    for(i in c(1:19,"X")){
      out2=scan1coef(qtl2apr[,as.character(i)],
                     expr[,j],
                     k[[i]],
                     sex)
      allmarker <- rbind(allmarker, out2$coef)
      print(paste( "coef scans on chr",i," of expression",j))
    }
    scanc[[j]]=allmarker
  }
  #formate coefficients into 3D array
  coef.array=abind(scanc[[1]], scanc[[2]], along=3)
  for(i in 3:length(scanc)){
    coef.array=abind(coef.array, scanc[[i]]$coef)
  }
  coef.array.t=aperm(coef.array, c(3,1,2))
  rownames(coef.array.t)=colnames(expr)
  
  #make the coef list
  coef=list()
  strains=data.frame(strain_id=probs$alleles, name=names(CCcolors))
  coef$strains=strains
  coef$coef=coef.array.t
  #insert coef into dataset
  dataset$coef=coef
}

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
  samples=data.frame(sample_id=rownames(genoprobs),
                     name=rownames(genoprobs), 
                     description=NA)
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
if(!is.null(expr)){
  dataset$expression$expression <- t(expr)
}


###save data##########
save(dataset, "eQTL_viewer.Rdata")



























