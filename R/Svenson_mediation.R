devtools::install_github("simecek/intermediate")
library(intermediate)

setwd("~/Projects/DO")
#load eQTL data and qlt2 data
load("data/DO478_ExprData_4eQTLMapping.Rdata")
load("data/Svenson_4eQTLMapping.Rdata")
#load phenotype scan1 data
load("results/7-26-16/Svenson_pheno_intdiet.RData")
load("results/7-26-16/Svenson_pheno_intsex.RData")
load("results/7-26-16/Svenson_pheno_add.RData")

#make interactive covariates
#for sex 
sexcovar=as.matrix(covar[,1])
colnames(sexcovar)=colnames(covar)[1]
#for diet
dietcovar=as.matrix(covar[,4])
colnames(dietcovar)=colnames(covar)[4]

#change names to the annotation data for mediation scan
names(annotations.rna.478)[4]="pos"

#do mediation scan with addcovar
med=vector("list")
for(i in 1:nrow(maxlod.add)){
  pp=maxlod.add$pheno[i]
  med[[pp]] <- mediation.scan(target = phenotype[,pp],
                             mediator = expr,
                             annotation = annotations.rna.478,
                             covar = covar,
                             qtl.geno = probs.478[,,maxlod.add$marker[which(
                               maxlod.add[,4]==pp)]])
}

#plot mediation scans and find genes that drop the lod score at the chromosome by 25%
#initialize a list
med_genes=vector("list")
#open pdf stream
pdf("results/7-21-16/mediation_x-20.pdf")
for(i in 1:nrow(maxlod.add)){
  name=maxlod.add$pheno[i]
  plot(med[[i]], main=paste(name,"| Chr", maxlod.add$chr[i]))
  media=subset(med[[i]], med[[i]]$Chr==maxlod.add$chr[i])
  print(i)
  if(nrow(media[which(media$LOD<(0.75*maxlod.add$LOD_scores[i])),])!=0){
    med_genes[[name]]=media[which(media$LOD<(0.75*maxlod.add$LOD_scores[i])),]
    print(med_genes[[name]])
  }
}
graphics.off()

