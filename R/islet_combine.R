.libPaths("/home/xchen/DO/scripts/library")
#if(!require(c("devtools", "RcppEigen"))) {
#  install.packages(c("devtools", "RcppEigen"),
#  repos = 'http://cran.rstudio.com/')
#  }
  
 library(devtools)
 install_github(paste0("rqtl/qtl2", c("geno",
  "scan", "plot", "convert")))
print("lib")
library(qtl2convert)
library(qtl2scan)
setwd("/home/xchen/DO/results/islet_DO/lod")
file=list.files()
length(file)
load(file[1], verbose=T)
scan1=out
load(file[2],verbose=T)
scan2=out
islet_eQTL_scan1=cbind.scan1(scan1, scan2)
i=1
for(f in file){
	load(f, verbose=T)
	islet_eQTL_scan1=cbind.scan1(islet_eQTL_scan1, out)
	print(i)
	i=i+1
	}
#print(head(islet_eQTL_scan1$lod))
print(dim(islet_eQTL_scan1$lod))
if(ncol(islet_eQTL_scan1$lod)>=47643){
for(f in file){
	system(paste("rm",f))
        print(paste("rm",f))
        }
}

print("saving file")
save(islet_eQTL_scan1, file="islet_eQTL_scan1.RData")

