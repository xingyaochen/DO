#Xingyao Chen + Clifton Jeffery
#7/05/16
#data formatting for heart DO data to do scans
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
print(args)
#path to your working directory
directory = args[1]
#path to the 4eQTL mapping data
dataPath = args[2]

#load in data
load(dataPath, verbose = T)

print("the printed names should correcspond to your args")

#names of the objects in your data
#name of phenotype data frame
phenotype = args[3]
#names of mRNA expression dataframe
expr = args[4]
#name of genotype probability dataframe
probs_doqtl = args[5]
#name of SNP dataframe
snps = args[6]
#path to save the final 4eQTLMapping data to
outName = args[7]

covar = cbind(phenotype$sex, phenotype$generation)
save(phenotype, snps, expr, probs_doqtl, covar, file = outName)
