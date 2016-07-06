#Xingyao Chen 
#6/30/16
#last edit: 7/06/16
#liver RNA expression data correlations with phenotypes
#linear models of liver RNA expression data against sex/diet covariates
setwd("~/Projects/DO")
#load the eQTL Rdata
load("data/DO478_ExprData_4eQTLMapping.Rdata")

#find and read the cleaned phenodata
file=list.files(path="data",pattern="clean")
pheno=read.csv(paste("data",file, sep="/"), header=T)[,-c(1:3)]

#subset so that only samples that have expr data are in pheno data
newpheno=pheno[match(names(probs.478[,1,1]), pheno$Sample),]
rownames(newpheno)=newpheno$Sample

#convert all phenotype data to numeric
datmat=matrix(as.numeric(as.matrix(newpheno)), nrow=nrow(newpheno))
colnames(datmat)=colnames(newpheno)
head(datmat)
rownames(datmat)=rownames(newpheno)

#load in some functions
source("scripts/myfunctions.R")

#find all correlations that are >0.3 
sigcor=data.frame()
for(i in 7:ncol(datmat)){
  cor=mycorr(expr.rna.478, datmat[,i])
  cor=cor[which(abs(cor)>0.3)]
  geneName=annotations.rna.478$Gene[match(names(cor), rownames(annotations.rna.478))]
  name=rep(colnames(datmat)[i], length(cor))
  sigcor=rbind(sigcor, cbind(geneName, name, cor))
  print(i)
  print(paste("length =" ,length(cor)))
}

# write.csv(sigcor, file="results/6-30-16/expr_cors.csv")



#find max cor for each phenotype
maxcor=data.frame()
for(i in 7:ncol(datmat)){
  cor=mycorr(expr.rna.478, datmat[,i])
  corMax=cor[which(cor==max(cor))]
  geneName=annotations.rna.478$Gene[match(names(corMax), rownames(annotations.rna.478))]
  pheno=rep(colnames(datmat)[i], length(corMax))
  maxcor=rbind(maxcor, cbind(geneName, pheno, corMax))
  print(i)
}
head(maxcor)
plot(as.factor(maxcor[,1]))

#export as csv
write.csv(maxcor, file="results/6-30-16/expr_maxcors.csv")


#find all correlations that are >90 percentile
percor=data.frame()
for(i in 7:ncol(datmat)){
  cor=mycorr(expr.rna.478, datmat[,i])
  ind=sort(abs(cor), decreasing = TRUE) #sort abs of cor in decreasing order
  top10=ind[1:ceiling(length(ind)/10)] #get the top 10% of correlation values
  cor=cor[abs(cor)%in%top10] #find in cors
  geneName=annotations.rna.478$Gene[match(names(cor), rownames(annotations.rna.478))] #use gene names rather than ensembl IDs
  name=rep(colnames(datmat)[i], length(cor)) #annotate cors
  percor=rbind(percor, cbind(geneName, name, cor))
  print(i)
  print(paste("length =" ,length(cor)))
}
#export as csv
write.csv(percor, file="results/7-06-16/expr_cors.csv")


#########################linear modeling#############################
#do lm and anovas for all gene expression values

#make some empty lists
mod.sexdiet=vector("list")
mod.sexXdiet=vector("list")
mod.sex=vector("list")
mod.diet=vector("list")
mod.null=vector("list")

#forloop through all 21454 expr RNAs
for(i in 1:ncol(expr.rna.478)){
  mod.sexdiet[[annotations.rna.478$Gene[i]]]=lm(expr.rna.478[,i]~newpheno$Sex+newpheno$Diet+newpheno$Gen, data=newpheno)
  mod.sexXdiet[[annotations.rna.478$Gene[i]]]=lm(expr.rna.478[,i]~Sex*Diet+Gen, data=newpheno)
  mod.sex[[annotations.rna.478$Gene[i]]]=lm(expr.rna.478[,i]~Sex+Gen, data=newpheno)
  mod.diet[[annotations.rna.478$Gene[i]]]=lm(expr.rna.478[,i]~Diet+Gen, data=newpheno)
  mod.null[[annotations.rna.478$Gene[i]]]=lm(expr.rna.478[,i]~Gen, data=newpheno)
  print(paste(i,"lm done!"))
}

#make some empty data frames
p.sex=data.frame()
p.diet=data.frame()
p.int=data.frame()

#forloop to bind all anova pvals into data frame
for(nm in names(mod.sex)){
  p.sex=rbind(p.sex, anova(mod.sexdiet[[nm]], mod.diet[[nm]])[2,6])
  p.diet=rbind(p.diet, anova(mod.sexdiet[[nm]], mod.sex[[nm]])[2,6])
  p.int=rbind(p.int, anova(mod.sexdiet[[nm]], mod.sexXdiet[[nm]])[2,6])
  print(paste(nm, "1 for-loop less to go!"))
}

#bind all into 1 giant dataframe
pval=data.frame(p.sex,p.diet, p.int)
rownames(pval)=names(mod.sex)
names(pval)=c('p.sex','p.diet','p.int')

#export as csv
write.csv(pval, file="../results/6-15-16/pval.csv")


#find all pvals that are <0.05 for significance
num=c()
for (i in 1:3)
  num=c(num, which(pval[,i]<0.05))
uni=unique(num)

#subset pval dataframe to get only the ones <0.05
pval.sig=pval[uni,]

#export pval.sig
write.csv(pval.sig, file="../results/6-15-16/pval_sig.csv")

#find lowest pval for each expr RNA
indx=vector()
min=vector()
for (i  in 1:nrow(pval.sig)){
  indx[i]=which.min(as.numeric(pval.sig[i,]))
  min[i]=min(as.numeric(pval.sig[i,]))
  
}
#lable as such
indx=gsub("1", "p.sex", x=indx)
indx=gsub("2", "p.diet", x=indx)
indx=gsub("3", "p.int", x=indx)

#make dataframe that has the genes, either sex, diet, or int, and the corresponding pval
rez=data.frame(genes=rownames(pval.sig), which=indx, pval=min)
#look as distribution
table(rez[,2])
#subset for just the p.int
interact=rez[which(rez[,2]=="p.int"),]

#export just the interact dataframe as csv
write.csv(interact, "results/6-30-16/expr_interactions.csv")











