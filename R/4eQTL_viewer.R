#Xingyao Chen + Clifton Jeffery
#7/05/16
#data formatting for inputting into eQTL viewer

.libPaths("/home/xchen/DO/scripts/library")
install.packages(c("devtools", "yaml",
                  "jsonlite", "data.table", "RcppEigen"), 
                 repos='http://cran.rstudio.com/')
library(devtools)
install_github(paste0("rqtl/qtl2", c("geno",
                                     "scan", "plot", "convert")))
library(qtl2geno)
library(qtl2scan)
library(qtl2plot)
library(qtl2convert)

load("../../../../hpcdata/gac/derived/Svenson_DO850/DO478_ExprData_4eQTLMapping.Rdata")
oldpheno=read.csv(paste("../data","DO_phenotype_data_clean.csv", sep="/"), header=T)

#eliminate the pheno data
# load("Churchill_project/Data/DO478_ExprData_4eQTLMapping.Rdata")
# oldpheno <- read.csv('Churchill_project/Data/data_restructured.csv',header=T)
newpheno <- oldpheno[match(names(probs.478[,1,1]), oldpheno$Sample),]
rownames(newpheno)=newpheno$Sample

###Variables##########
nameOfExperiment <- "DO478"
phenotypeDescriptions <- read.csv("../data/svenson850_phenotype_descr.csv")
phenotype <- newpheno
genotype <- NULL
expression <- expr.rna.478
geneAnnotations <- annotations.rna.478
snps <- snps.64K
probs <- probs.478

######################

phenotype$urine.microalbumin1 <- as.numeric(as.character(phenotype$urine.microalbumin1))


###features###########
if(!is.null(geneAnnotations)){
  dataset=list()
  features=list()
  #annotations are the gene annotatons with ensemblID, chr, location, gene description,
  #and gene name
  features$feature_id=geneAnnotations$EnsemblID
  features$group_id=geneAnnotations$EnsemblID
  features$chrom=geneAnnotations$Chr
  features$locations=geneAnnotations$Start.Mbp
  features$name=geneAnnotations$Gene
  features$description=geneAnnotations$Gene.Biotype
  dataset$features=features
  
  
  ###markers###########
  if(!is.null(snps)){
    #snps are set as the  snp info
    markers=list()
    markers$marker_id=snps$SNP_ID
    markers$chrom=snps$Chr
    markers$location=snps$Mb_NCBI38
    #put everything into dataset
    
    dataset$markers=markers
  }
  
  ###QTL scans##########
  if(!is.null(probs) && !is.null(snps) && !is.null(expression)){
    if (class(probs) != "calc_genoprob"){
      probs=probs_doqtl_to_qtl2(probs, map=snps,
                                chr_column = "Chr", pos_column = "cM",
                                marker_column = "SNP_ID")
    }
    qtl2apr = genoprob_to_alleleprob(probs)
    
    #get kinship and sex covars
    k <- calc_kinship(qtl2apr, "loco", cores=16)
    sex <- (newpheno$Sex=="M")*1
    names(sex)=rownames(newpheno)
    print("scan1-ing")
    lod=list()
    lodvalue=data.frame(nrow=nrow(snps))
    for(i in 1:ncol(expression)){
      pheno=as.matrix(expression[,i])
      colnames(pheno)=colnames(expression)[i]
      out <- scan1(qtl2apr, pheno, k, sex, cores=16)
      print(paste("finished scan1 expr", i))
      lodvalue=cbind(lodvalue, out$lod)
    }
    
    
    lod$lod=lodvalue[,-1]
    dataset$lod=lod
  }
  
  
  ###Coef Scans#########
  if(!is.null(probs) && !is.null(snps) && !is.null(expression)){
    #do the coef thing
    scanc=vector("list")
    for(j in 1:ncol(expression)){
      allmarker=data.frame()
      for(i in c(1:19,"X")){
        out2=scan1coef(qtl2apr[,as.character(i)],
                       expression[,j],
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
    rownames(coef.array.t)=colnames(expression)
    
    #make the coef list
    coef=list()
    coef$name=names(CCcolors)
    coef$strains=colnames(probs)
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
    samples$sample_id=rownames(phenotype)
    samples$name=rownames(phenotype)
    dataset$samples = samples
  }
  
  
  ###phenotypes#########
  if(!is.null(phenotypeDescriptions) && !is.null(phenotype)){
    factors <- data.frame()
    for (i in 1:length(phenotype)) {
      index <- grep(colnames(phenotype)[i], phenotypeDescriptions$name)
      insert <- cbind(as.character(phenotypeDescriptions$name[index]),
                      as.character(phenotypeDescriptions$name[index]),
                      as.character(phenotypeDescriptions$description[index]))
      factors <- rbind(factors, insert)
    }
    names(factors) <- c("factor_id", "name", "description")
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
  save(dataset, "../data/viewer.Rdata")
  
  
  
  
  
  
  
  
    
    
    
    
  
  
  
  
  
  