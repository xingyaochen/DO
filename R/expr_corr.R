setwd("~/Projects/DO")
load("data/DO478_ExprData_4eQTLMapping.Rdata")
file=list.files(path="data",pattern="clean")
pheno=read.csv(paste("data",file, sep="/"), header=T)[,-c(1:3)]

newpheno=pheno[match(names(probs.478[,1,1]), pheno$Sample),]
rownames(newpheno)=newpheno$Sample

datmat=matrix(as.numeric(as.matrix(newpheno)), nrow=nrow(newpheno))
colnames(datmat)=colnames(newpheno)
head(datmat)
rownames(datmat)=rownames(newpheno)

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

write.csv(maxcor, file="results/6-30-16/expr_maxcors.csv")



mod.sexdiet=vector("list")
mod.sexXdiet=vector("list")
mod.sex=vector("list")
mod.diet=vector("list")
mod.null=vector("list")

for(i in 1:ncol(expr.rna.478)){
  mod.sexdiet[[annotations.rna.478$Gene[i]]]=lm(expr.rna.478[,i]~newpheno$Sex+newpheno$Diet+newpheno$Gen, data=newpheno)
  mod.sexXdiet[[annotations.rna.478$Gene[i]]]=lm(expr.rna.478[,i]~Sex*Diet+Gen, data=newpheno)
  mod.sex[[annotations.rna.478$Gene[i]]]=lm(expr.rna.478[,i]~Sex+Gen, data=newpheno)
  mod.diet[[annotations.rna.478$Gene[i]]]=lm(expr.rna.478[,i]~Diet+Gen, data=newpheno)
  mod.null[[annotations.rna.478$Gene[i]]]=lm(expr.rna.478[,i]~Gen, data=newpheno)
  print(paste(i,"lm done!"))
}

p.sex=data.frame()
p.diet=data.frame()
p.int=data.frame()
for(nm in names(mod.sex)){
  p.sex=rbind(p.sex, anova(mod.sexdiet[[nm]], mod.diet[[nm]])[2,6])
  p.diet=rbind(p.diet, anova(mod.sexdiet[[nm]], mod.sex[[nm]])[2,6])
  p.int=rbind(p.int, anova(mod.sexdiet[[nm]], mod.sexXdiet[[nm]])[2,6])
  print(paste(nm, "1 for-loop less to go!"))
}

pval=data.frame(p.sex,p.diet, p.int)
rownames(pval)=names(mod.sex)
names(pval)=c('p.sex','p.diet','p.int')
setwd("..")
getwd()
write.csv(pval, file="results/6-15-16/pval.csv")
num=c()
for (i in 1:3)
  num=c(num, which(pval[,i]<0.05))

uni=unique(num)
pval.sig=pval[uni,]
#write.csv(pval.sig, file="results/6-15-16/pval_sig.csv")
head(pval.sig)
indx=vector()
min=vector()
for (i  in 1:nrow(pval.sig)){
  indx[i]=which.min(as.numeric(pval.sig[i,]))
  min[i]=min(as.numeric(pval.sig[i,]))
  
}
indx=gsub("1", "p.sex", x=indx)
indx=gsub("2", "p.diet", x=indx)
indx=gsub("3", "p.int", x=indx)

rez=data.frame(genes=rownames(pval.sig), which=indx, pval=min)
table(rez[,2])
interact=rez[which(rez[,2]=="p.int"),]

write.csv(interact, "results/6-30-16/expr_interactions.csv")
dim(na.omit(expr.rna.478))
pc_expr=princomp(x=expr.rna.478, cor=T)


