#Xingyao Chen
#6/23/16
#PCA for all phenotypes
setwd("~/Projects/DO/data")
file=list.files(pattern="clean")
data=read.csv(file, header=T)

data=newpheno
source("../scripts/myfunctions.R")
#reomve all NAs from dataset
naList <- data.frame()
for(i in (grep("Coat.Color", colnames(data))+1):dim(data)[2]){
  insert <- cbind(colnames(data[i]), length(data[,i]) - length(na.omit(data[,i])))
  naList <- rbind(naList, insert)
}

ls <- subset(naList, as.numeric(as.character(naList[,2])) > 100)
index=as.numeric(match(as.character(ls[,1]),as.character(names(data))))
data.reduced <- data[,-index]
data.reduced <- na.omit(data.reduced)
head(data.reduced)

#change covariates back to characters
data.reduced[which(data.reduced[,4]==1),4]="F"
data.reduced[which(data.reduced[,4]==0),4]="M"
data.reduced[which(data.reduced[,7]==1),7]="hf"
data.reduced[which(data.reduced[,7]==0),7]="chow"

#rz transform data
for(i in c(9:ncol(data.reduced))){
  data.reduced[,i]=rz.transform(as.numeric(data.reduced[,i]))
}


#look
head(data.reduced)
#get the grouping, group by sex and diet together
data.reduced.sexDiet=paste(data.reduced$Sex, data.reduced$Diet, sep="_")

#princomp
pc <- princomp(x = data.reduced[,9:ncol(data.reduced)], cor = TRUE)
predict(pc, 
        newdata=tail(data.reduced, 2))
#find phenotypes for comp 1-5
pc.phenos=vector("list")
for(i in 1:5){
  pc.phenos[[i]]=pc$loadings[which(abs(pc$loadings[,i])>0.1),i]
}

#plot them out as biplots
ggbiplot(pc, choices=c(1,2),groups=data.reduced.sexDiet, 
         ellipse = TRUE,alpha=0.2,
         circle = TRUE, varname.size = 3,
         varname.abbrev = TRUE, size=0.5)+
  ggtitle("RZ transformed Phenotypes")+
  theme_classic()

ggbiplot(pc, choices=c(1,3),groups=data.reduced.sexDiet, 
              ellipse = TRUE,alpha=0.2,
              circle = TRUE, varname.size = 3,
              varname.abbrev = TRUE, size=0.5)+
  ggtitle("RZ transformed Phenotypes")+
  theme_classic()

ggbiplot(pc, choices=c(2,3),groups=data.reduced.sexDiet, 
         ellipse = TRUE,alpha=0.2,
         circle = TRUE, varname.size = 3,
         varname.abbrev = TRUE, size=0.5)+
  ggtitle("RZ transformed Phenotypes")+
  theme_classic()

ggbiplot(pc, choices=c(3,4),groups=data.reduced.sexDiet, 
         ellipse = TRUE,alpha=0.2,
         circle = TRUE, varname.size = 3,
         varname.abbrev = TRUE, size=0.5)+
  theme_classic()



ggbiplot(pc, choices=c(2,1),groups=data.reduced.sexDiet, 
         ellipse = TRUE,alpha=0.2,
         circle = TRUE, varname.size = 3,
         varname.abbrev = TRUE, size=0.5)+
  theme_classic()


#plot histogram of pc1 scores, faceted by sex and diet
pc.data=cbind(data.reduced[,c(4,7)], as.data.frame(pc$scores))
pc.data$Sex

ggplot(data=pc.data, 
       aes(pc.data[,3], fill=c(Sex)))+
  geom_histogram(bins=60, alpha=.5, position="identity")+
  facet_grid(Diet~Sex)+
  labs(x="Comp.1")

#make linear models
mod.sexdiet=vector("list")
mod.sexXdiet=vector("list")
mod.sex=vector("list")
mod.diet=vector("list")
mod.null=vector("list")

for(i in 1:ncol(pc$scores)){
  mod.sexdiet[[i]]=lm(pc$scores[,i]~data.reduced$Sex+data.reduced$Diet+data.reduced$Gen, data=data.reduced)
  mod.sexXdiet[[i]]=lm(pc$scores[,i]~Sex*Diet+Sex+Diet+Gen, data=data.reduced)
  mod.sex[[i]]=lm(pc$scores[,i]~data.reduced$Sex+data.reduced$Gen, data=data.reduced)
  mod.diet[[i]]=lm(pc$scores[,i]~data.reduced$Diet+data.reduced$Gen, data=data.reduced)
  mod.null[[i]]=lm(pc$scores[,i]~data.reduced$Gen, data=data.reduced)
}

p.sex=data.frame()
p.diet=data.frame()
p.int=data.frame()
for(nm in 1:length(mod.sex)){
  p.sex=rbind(p.sex, anova(mod.sexdiet[[nm]], mod.diet[[nm]])[2,6])
  p.diet=rbind(p.diet, anova(mod.sexdiet[[nm]], mod.sex[[nm]])[2,6])
  p.int=rbind(p.int, anova(mod.sexdiet[[nm]], mod.sexXdiet[[nm]])[2,6])
}

pval=data.frame(p.sex,p.diet, p.int)
rownames(pval)=colnames(pc$scores)
names(pval)=c("p.sex","p.diet","p.int")


num=c()
for (i in 1:3)
  num=c(num, which(pval[,i]<0.05))

uni=sort(unique(num))
pval.sig=pval[uni,]
dim(pval.sig)

write.csv(pval.sig, file="../results/6-23-16/pca_sig_pvals.csv")

