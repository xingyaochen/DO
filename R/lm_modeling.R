#Xingyao Chen
#6/15/16 
#phenotype basic lm modeling
#
setwd("~/Projects/DO/data")
file=list.files()
#read phenotype data
data=read.csv(file, header=T)

#make empty lists
mod.sexdiet=vector("list")
mod.sexXdiet=vector("list")
mod.sex=vector("list")
mod.diet=vector("list")
mod.null=vector("list")

#do linear regression with phenotype against covariates
for(i in 7:ncol(data)){
  mod.sexdiet[[names(data)[i]]]=lm(data[,i]~data$Sex+data$Diet+data$Gen, data=data)
  mod.sexXdiet[[names(data)[i]]]=lm(data[,i]~Sex*Diet+Gen, data=data)
  mod.sex[[names(data)[i]]]=lm(data[,i]~data$Sex+data$Gen, data=data)
  mod.diet[[names(data)[i]]]=lm(data[,i]~data$Diet+data$Gen, data=data)
  mod.null[[names(data)[i]]]=lm(data[,i]~data$Gen, data=data)
}

#make empty data frames
p.sex=data.frame()
p.diet=data.frame()
p.int=data.frame()

#add anova pvals into data frames
for(nm in names(mod.sex)){
  p.sex=rbind(p.sex, anova(mod.sexdiet[[nm]], mod.diet[[nm]])[2,6])
  p.diet=rbind(p.diet, anova(mod.sexdiet[[nm]], mod.sex[[nm]])[2,6])
  p.int=rbind(p.int, anova(mod.sexdiet[[nm]], mod.sexXdiet[[nm]])[2,6])
}
#combine into 1 giant data frame
pval=data.frame(p.sex,p.diet, p.int)
rownames(pval)=names(mod.sex)
names(pval)=c('p.sex','p.diet','p.int')

getwd()
#export as csv
write.csv(pval, file="../results/6-15-16/pval.csv")

#find significant pvals
num=c()
for (i in 1:3)
  num=c(num, which(pval[,i]<0.05))

uni=unique(num)

#subset pval dataframe to get only significant ones as export as csv
pval.sig=pval[uni,]
write.csv(pval.sig, file="results/6-15-16/pval_sig.csv")


#find the mininum pval
indx=vector()
min=vector()
for (i  in 1:nrow(pval.sig)){
  indx[i]=which.min(as.numeric(pval.sig[i,]))
  min[i]=min(as.numeric(pval.sig[i,]))
  
}
indx=gsub("1", "p.sex", x=indx)
indx=gsub("2", "p.diet", x=indx)
indx=gsub("3", "p.int", x=indx)

#make a dataframe with phenotype, which covar, and the p-value
rez=data.frame(phenotypes=rownames(pval.sig), which=indx, pval=min)
dim(rez)
#export as csv
write.csv(rez, file="../results/6-15-16/pval_min.csv")

min.diet=rez[grep("*diet", rez[,2]),]
min.sex=rez[grep("*sex", rez[,2]),]
# write.csv(min.diet, file="results/6-15-16/pval_min_diet.csv")
# write.csv(min.sex, file="results/6-15-16/pval_min_sex.csv")
# 













