.libPaths("/home/xchen/DO/scripts/library")
install.packages("qtlcharts",  repos = 'http://cran.rstudio.com/')
library(qtlcharts)
getwd()
setwd("/home/xchen/DO")
load("data/Svenson_4eQTLMapping.RData")

corpheno=iplotCorr(phenotype[,7:ncol(phenotype)],covar[,'sex'], reorder=TRUE,
          chartOpts=list(cortitle="Correlation matrix",
                         scattitle="Scatterplot | Sex"))

save(corpheno, file="results/corpheno.RData")

corexpr=iplotCorr(expr,covar[,'sex'], reorder=TRUE,
          chartOpts=list(cortitle="Correlation matrix",
                         scattitle="Expression Scatterplot | Sex"))
save(corexpr, file="results/corexpr.RData")

