#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
directory=args[1]
outFile=args[2]
.libPaths(paste(directory, "scripts/library", sep="/"))
library(qtl2scan)
setwd(directory)
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
	all_eQTL_scan1=cbind.scan1(islet_eQTL_scan1, out)
	print(i)
	i=i+1
	}
#print(head(islet_eQTL_scan1$lod))
print(dim(all_eQTL_scan1$lod))
for(f in file){
	system(paste("rm",f))
        print(paste("rm",f))
        }

print("saving file")
save(all_eQTL_scan1, file=outFile)

