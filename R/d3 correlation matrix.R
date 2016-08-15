install.packages("qtlcharts")
library(qtlcharts)
setwd("~/Projects/DO/data")
file=list.files(pattern="clean")
data=read.csv(file, header=T)
head(data)
dim(data)
#make the phenotypes numeric
pheno=matrix(as.numeric(as.matrix(data[,-c(1:9)])), nrow=nrow(data[,-c(1:9)]))
colnames(pheno)=names(data[,-c(1:9)])
rownames(pheno)=data$Sample
head(pheno)

#factor sex and diet 
sexDiet=paste(data$Sex, data$Diet, sep="_")
sexDiet=factor(sexDiet)

#plot the correlation matrix factored by sex and diet
iplotCorr(pheno,sexDiet, reorder=TRUE,
          chartOpts=list(cortitle="Correlation matrix",
                         scattitle="Scatterplot | Sex and Diet",
                         scatcolors=c("#FA5858","#86B404","#01DFD7","#BE81F7")
                         ))

iplotCorr(pheno,data$Sex, reorder=TRUE,
          chartOpts=list(cortitle="Correlation matrix",
                         scattitle="Scatterplot | Sex",
                         scatcolors=c("red","#00BFFF")
          ))

iplotCorr(pheno,data$Diet, reorder=TRUE,
          chartOpts=list(cortitle="Correlation matrix",
                         scattitle="Scatterplot | Diet",
                         scatcolors=c("#298A08","#FF8000")
          ))


y =x
###plot the legends because idk ow to do them on the iplot
#sexdiet legend
plot(0, 0, type = "n", yaxt="n", xaxt="n", xlab="", ylab="")
legend(-1,1, as.character(unique(sexDiet)), pch = 21,
       pt.bg = c("#FA5858","#86B404","#01DFD7","#BE81F7"),
       bty = "n",
       y.intersp	=2)

#sex legend
plot(0, 0, type = "n", yaxt="n", xaxt="n", xlab="", ylab="")
legend(-1,1, as.character(unique(data$Sex)), pch = 21,
       pt.bg =c("red","#00BFFF"),
       inset = 1,
       bty = "n",
       y.intersp	=2)

#diet legend
plot(0, 0, type = "n", yaxt="n", xaxt="n", xlab="", ylab="")
legend(-1,1, as.character(unique(data$Diet)), pch = 21,
       pt.bg =c("#298A08","#FF8000"),
       inset = 1,
       bty = "n",
       y.intersp	=2)