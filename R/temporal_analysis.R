#Xingyao Chen
#6/17/16
#Plotting phenotypes vs time
library(ggplot2)
library(corrplot)
library(plyr)
library(reshape2)
setwd("~/Projects/DO/data")
file=list.files(pattern="svenson.phenotypes.final")
data=read.csv(file, header=T)
dim(data)
#scale transform data
# for(i in c(7:150)){
#   data[,i]=scale(as.numeric(data[,i]))
# }

setwd("..")

#find all the phenotypes that have 1 and 2 at end, those are temporal phenotypes of interest
indtemp=c(grep("*1",names(data), value=F),grep("*2",names(data), value=F))
subset=data[,c(3:5,sort(indtemp))]
#take out the body weight data, will deal with those later
minus=grep("BW.*", names(subset), value=F)
timDat=subset[,-c(minus)]

#these are the phenotypes with 2 time points
names(timDat)

#reformat data and delet the 1 and 2
meltDat=melt(timDat, id.vars=c("Sex","Gen","Diet"))
#delete the 1 and 2 at end to make phenotypes same
#remove all [.] first
meltDat$variable=gsub("[.]","", meltDat$variable)
#gsub back in ".1" and ".2"
meltDat$variable=gsub("2",".2", meltDat$variable)
meltDat$variable=gsub("1",".1", meltDat$variable)
#split by [.], then unsplit and add back into frame
ls=unlist(strsplit(meltDat$variable,"[.]"))
ls.mat=matrix(ls, nrow=2)
dim(ls.mat)
meltDat$variable=ls.mat[1,]
meltDat$week=ls.mat[2,]

#split data into a list of 49 phenotypes
meltList=split(meltDat, meltDat$variable)
tail(meltList[[1]])
#remove all the phenotyeps w/o 2 discreet time points (again)
list=c()
for(i in 1:49){
  if(length(unique(meltList[[i]]$week))!=2){
    list=c(list, i)
  }
}
meltList=meltList[-list]
#new length
length(meltList)

setwd("~/Projects/DO")
#open a pdf stream
pdf("results/6-20-16/temporal_plots_lines_transformed.pdf")
for(i in 1:length(meltList)){
  #define the weeks
  meltList[[i]]$value=scale(as.numeric(meltList[[i]]$value))
  a=desc[grep(pattern=paste(meltList[[i]]$variable, "*", sep=""), desc[,1]),][1,9]
  b=desc[grep(pattern=paste(meltList[[i]]$variable, "*", sep=""), desc[,1]),][2,9]
  #plot phenotype by diet
p1=ggplot(na.omit(meltList[[i]]), aes(y=as.numeric(value), x=week, group=Diet, col=Diet))+
  geom_smooth()+
  geom_jitter(width=0.08, size=0.1)+
  theme_classic()+
  labs(y="value")+
  scale_x_discrete(labels=c(a,b))+
  theme(legend.title=element_blank(), 
        plot.title=element_text(size=9),
        legend.position="top") +
  ggtitle(paste(names(meltList)[i],"| Diet"))

#plot phenotype by Sex
p2=ggplot(na.omit(meltList[[i]]), aes(y=as.numeric(value), x=week, group=Sex, col=Sex))+
  geom_smooth()+
  geom_jitter(width=0.08, size=0.1)+
  theme_classic()+
  labs(y="value")+
  scale_x_discrete(labels=c(a,b))+
  theme(legend.title=element_blank(), 
        plot.title=element_text(size=9),
        legend.position="top") +
  ggtitle(paste(names(meltList)[i],"| Sex"))

#plot phenotype by Gen (took away)
# p1=ggplot(na.omit(meltList[[i]]), aes(y=as.numeric(value), x=week, group=Gen, col=Gen))+
#   geom_smooth()+
#   theme_classic()+
#   labs(y="value")+
#   theme(legend.title=element_blank(), 
#         axis.ticks = element_blank(), 
#         axis.text.y = element_blank(),
#         plot.title=element_text(size=9),
#         legend.position="top") +
#   ggtitle(paste(unique(meltDat$variable)[i],"| Gen"))

print(multiplot(p1, p2, cols=2))
}

#close PDF
dev.off()


coef.diet=data.frame()
coef.sex=data.frame()

for(i in 1:length(meltList)){
  #do lm with each phenotype against week, subset by hf
  coef.hf=lm(as.numeric(value)~week, data=na.omit(meltList[[i]]), Diet=="hf")$coefficients[2]
  #subset by chow
  coef.chow=lm(as.numeric(value)~week, data=na.omit(meltList[[i]]), Diet=="chow")$coefficients[2]
  #subset by Females
  coef.f=lm(as.numeric(value)~week, data=na.omit(meltList[[i]]), Sex=="F")$coefficients[2]
  #subset by males
  coef.m=lm(as.numeric(value)~week, data=na.omit(meltList[[i]]), Sex=="M")$coefficients[2]
  #make into dataframe
  this=data.frame(chow=coef.chow, hf=coef.hf)
  coef.diet=rbind(coef.diet, this)
  that=data.frame(Female=coef.f, Male=coef.m)
  coef.sex=rbind(coef.sex, that)
}
#add back phenotype names
rownames(coef.diet)=names(meltList)
rownames(coef.sex)=names(meltList)
#export to csv
write.csv(coef.diet, file="results/6-20-16/coef_sub_by_diet_transformed.csv")
write.csv(coef.sex, file="results/6-20-16/coef_sub_by_sex_transformed.csv")

#find all coefs that have opposite signs between m/f or hf/chow
putout.diet=c()
putout.sex=c()
for(i in 1:nrow(coef.diet)){
  #extract all coef that have a difference of 0.1 and have opposite signs (+/-)
  if(abs(coef.diet[i,1]-coef.diet[i,2])>=0.05){
    if(coef.diet[i,1]>0 & coef.diet[i,2]<0 | coef.diet[i,1]<0 & coef.diet[i,2]>0){
      putout.diet=c(putout.diet, i)
  }}
  if(abs(coef.sex[i,1]-coef.sex[i,2])>=0.05){
    if(coef.sex[i,1]>0 & coef.sex[i,2]<0 | coef.sex[i,1]>0 & coef.sex[i,2]<0){
      putout.sex=c(putout.sex, i)
  }
  }
}
#these
coef.diet[putout.diet,]
#         chow         hf
# CHOL -0.7524335   3.018110
# GLDH -2.3155801   4.489686
# PLT   0.8714324 -15.905339

#transformed data
#           chow          hf
# CHOL -0.023685326  0.09500496
# GLDH -0.210269739  0.40769268
# PLT   0.002625295 -0.04791675

coef.sex[putout.sex,]
#                     Female        Male
# CHOL              0.94011876 -0.51639288
# HCT               0.62590430 -0.13796652
# urinecreatinine   2.99890513 -1.19205594
# urinemicroalbumin 0.09497205 -0.03962201

#transformed data
                      #Female        Male
# ACR               0.1199830 -0.06370763
# HCT               0.1441325 -0.03177077
# urinecreatinine   0.1174775 -0.04669698
# urinemicroalbumin 0.1094498 -0.04566209



#now deal with body weight data
#grep all the BW data
data[which(data[,4]==1),4]="F"
data[which(data[,4]==0),4]="M"
data[which(data[,7]==1),7]="hf"
data[which(data[,7]==0),7]="chow"

indBW=grep("BW.*",names(data))
BW.dat=data[,c(3:7, indBW)]
#look
head(BW.dat)
dim(BW.dat)

#reformat and take out NAs
BWmelt=melt(BW.dat, id.vars = c("Sample","Sex","Gen","parity", "Diet"))
BWmelt.clean=na.omit(BWmelt)
#dcast(BWmelt.clean, Sex~variable, fun.aggregate = median)


#plot phenotype by diet
#plot phenotype by diet
setwd("/Users/xchen/Projects/DO/results")
p4=ggplot(BWmelt, aes(x=as.numeric(variable), y=as.numeric(value), col=Sex, group=Sex))+
  geom_smooth()+
  geom_jitter(width=0.5, size=0.01)+
  ggtitle("Body weight | Sex")+
  theme_classic()+
  labs(y="grams", x="Age")
p5=ggplot(BWmelt, aes(x=as.numeric(variable), y=as.numeric(value), col=Diet, group=(Diet)))+
  geom_smooth()+
  geom_jitter(width=0.5, size=0.01)+
  ggtitle("Body weight | Diet")+
  theme_classic()+
  labs(y="grams", x="Age")
p6=ggplot(BWmelt, aes(x=as.numeric(variable), y=as.numeric(value), col=Gen, group=Gen))+
  geom_smooth()+
  geom_jitter(width=0.5, size=0.01)+
  ggtitle("Body weight | Gen")+
  theme_classic()+
  labs(y="grams", x="Age")

p7=ggplot(BWmelt, aes(x=as.numeric(variable), y=as.numeric(value), col=parity, group=parity))+
  geom_smooth()+
  geom_jitter(width=0.5, size=0.01)+
  ggtitle("Body weight | Parity")+
  theme_classic()+
  labs(y="grams", x="Age")

pdf("6-20-16/BW_over_time.pdf")
print(p4)
print(p5)
print(p6)
print(p7)
dev.off()
graphics.off()
