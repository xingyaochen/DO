#Xingyao Chen + Clifton Jeffery
#7/05/16
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

#specify the location to which libraries are installed (make a folder named library in 
#your working directory)

.libPaths(paste(directory,"library", sep="/"))

#install packages if not already installed

if(!require(devtools))
  install.packages(c("devtools", "yaml",
                     "jsonlite", "data.table", "RcppEigen"), 
                   repos='http://cran.rstudio.com/')
library(devtools)


if(!require("biomaRt")){
  source("https://bioconductor.org/biocLite.R")
  biocLite("biomaRt")
}


if(!require("qtl2geno")|!require("qtl2scan")|!require("qtl2plot")|
   !require("qtl2convert"))
install_github(paste0("rqtl/qtl2", c("geno",
                                     "scan", "plot", "convert")))

library(biomaRt)
library(qtl2scan)
library(qtl2plot)
library(qtl2convert)

#load in the RData for filepath
load(filePath)
#set current working directory
setwd(directory)

#check to see if sample numbers match
stopifnot(rownames(expr.mrna)==rownames(pheno)|
         rownames(expr.mrna)==rownames(probs))


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
    k <- calc_kinship(qtl2apr, "loco", cores=4)
    sex <- (pheno$sex=="M")*1
    names(sex)=rownames(pheno)
    print("scan1-ing")
    
    pheno=as.matrix(expr[,1])
    rownames(pheno)=covar[,1]
    colnames(pheno)=colnames(expr)[1]
    allscan1=scan1(qtl2apr, pheno, k, sex, cores=4)
    
    lod=list()
    lodvalue=data.frame(nrow=nrow(snps))
    for(i in num:c(num+999)){
      pheno=as.matrix(expr[,i])
      rownames(pheno)=covar[,1]
      colnames(pheno)=colnames(expr)[i]
      out <- scan1(qtl2apr, pheno, k, sex, cores=4)
      print(paste("finished scan1 expr", i))
      lodvalue=cbind(lodvalue, out$lod)
      write.csv(lodvalue, paste("lod", num, ".csv", sep=""))
      if(i>1){
        allscan1=cbind.scan1(allscan1, out)
        save(allscan1, "allscan1.RData")
      }  
      if(num==ncol(expr)){
        break
      }
    }
    lod$lod=lodvalue[,-1]
    dataset$lod=lod
  }

if(num==floor(ncol(expr)/1000)*1000){
  
    
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
  if(!is.null(expression)){
    dataset$expression$expression <- t(expression)
  }

  ###save data##########
  save(dataset, "viewer.Rdata")
  
  
  
  
  
  
  
  
    
    
    
    
  
  
  
  
  
  