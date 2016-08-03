load(("data/DO478_ExprData_4eQTLMapping.Rdata"))
load("data/Svenson_4eQTLMapping.Rdata")
annot=annotations.rna.478
probs_doqtl=probs.478
phenotype_descr=read.csv("data/svenson850_phenotype_descr.csv")
save(covar, expr, phenotype_descr, phenotype, snps, k, probs_doqtl, probs, 
     file="data/Svenson_4eQTLMapping.Rdata")
