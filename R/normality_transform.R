#Xingyao Chen
#6/22/16
#transforming pheno data to see which scale works best
setwd("~/Projects/DO/data")
#read in pheo data
file = list.files(pattern = "final")
data = read.csv(file, header = T)

source("../scripts/my_funtions.R")
source("http://peterhaschke.com/Code/multiplot.R")


pheno = data[, c(1:6)]
mydata = data[, -c(1:6)]
#RankZ transform data
data.rz = data
for (i in c(7:150)) {
  data.rz[, i] = rz.transform(as.numeric(data[, i]))
}
#name
names(data.rz) = names(data)
#
#Log transform data
data.log = data

for (i in c(7:150)) {
  data.log[, i] = log(as.numeric(data[, i]))
}
dim(data.log)


#######################
#plot all pheno histograms using non transformed, normal data
########################
#
pdf("../results/6-22-16/pheno_histogram.pdf")
for (i in 7:150) {
  p1 = ggplot(data, aes(x = as.numeric(data[, i]), fill = Sex)) +
    geom_histogram(bins = 60,
                   alpha = .5,
                   position = "identity") +
    facet_grid(Diet ~ .) +
    theme(legend.position = 'top') +
    labs(x = "") +
    ggtitle(names(data)[i])
  
  
  p2 = ggplot(data, aes(x = as.numeric(data[, i]), fill = Diet)) +
    geom_histogram(bins = 60,
                   alpha = .5,
                   position = "identity") +
    scale_fill_brewer(palette = "Set1") +
    facet_grid(Sex ~ .) +
    theme(legend.position = 'top') +
    labs(x = "") +
    ggtitle(names(data)[i])
  
  print(multiplot(p1, p2, cols = 2))
}
graphics.off()


#######################
#plot all pheno histograms using log transformed data
########################

#save all plots to pdf
pdf("../results/6-22-16/pheno_log_histogram_.pdf")
for (i in 7:150) {
  #plot histogram, colored by Sex, faceted by Diet
  p1 = ggplot(data.log, aes(x = as.numeric(data.log[, i]), fill = Sex)) +
    geom_histogram(bins = 60,
                   alpha = .5,
                   position = "identity") +
    facet_grid(Diet ~ .) +
    labs(x = "") +
    theme(legend.position = 'top') +
    ggtitle(names(data.log)[i])
  
  #plot histogram, colored by Diet, faceted by Sex
  p2 = ggplot(data.log, aes(x = as.numeric(data.log[, i]), fill = Diet)) +
    geom_histogram(bins = 60,
                   alpha = .5,
                   position = "identity") +
    scale_fill_brewer(palette = "Set1") +
    facet_grid(Sex ~ .) +
    labs(x = "") +
    theme(legend.position = 'top') +
    ggtitle(names(data.log)[i])
  
  print(multiplot(p1, p2, cols = 2))
}
graphics.off()


#######################
#plot all pheno histograms using rz transformed data
########################
pdf("results/6-22-16/pheno_rz_histogram_.pdf")
for (i in 7:150) {
  p1 = ggplot(data.rz, aes(x = as.numeric(data.rz[, i]), fill = Sex)) +
    geom_histogram(bins = 60,
                   alpha = .5,
                   position = "identity") +
    facet_grid(Diet ~ .) +
    labs(x = "") +
    theme(legend.position = 'top') +
    ggtitle(names(data.rz)[i])
  
  
  p2 = ggplot(data.rz, aes(x = as.numeric(data.rz[, i]), fill = Diet)) +
    geom_histogram(bins = 60,
                   alpha = .5,
                   position = "identity") +
    scale_fill_brewer(palette = "Set1") +
    facet_grid(Sex ~ .) +
    labs(x = "") +
    theme(legend.position = 'top') +
    ggtitle(names(data.rz)[i])
  
  print(multiplot(p1, p2, cols = 2))
}
graphics.off()


#do lillie normality tests, pvals <0.1 is considered significat
ntest = data.frame()
ntest.rz = data.frame()
ntest.log = data.frame()
for (i in 7:ncol(data)) {
  x = lillie.test(as.numeric(data[, i]))$p.value
  rz = lillie.test(as.numeric(data.rz[, i]))$p.value
  log = sf.test(data.log[, i])$p.value
  ntest = rbind(ntest, x)
  ntest.rz = rbind(ntest.rz, rz)
  ntest.log = rbind(ntest.log, log)
}
names(ntest) = "pval"
names(ntest.log) = "pval"
names(ntest.rz) = "pval"
ntest$phenotypes = names(data)[7:150]
ntest.rz$phenotypes = names(data)[7:150]
ntest.log$phenotypes = names(data)[7:150]
head(ntest)

ntest.sig = ntest[which(ntest[, 1] < 0.05), ]
ntest.rz.sig = ntest[which(ntest.rz[, 1] < 0.05), ]
ntest.log.sig = nest[which(ntest.log[, 1] < 0.05), ]
nrow(ntest.sig)
nrow(ntest.rz.sig)
nrow(ntest.log.sig)

write.csv(ntest, file = "normality-test.csv")
write.csv(ntest.rz, file = "normality-test_rz.csv")
write.csv(ntest.log, file = "normality-test_log.csv")

pdf("results/6-22-16/qqnorm_compare.pdf")
par(mfrow = c(1, 3))
for (i in 7:150) {
  qqnorm(as.numeric(data[, i]), main = paste("Normal Q-Q Plot",
                                             names(data)[i]))
  qqline(as.numeric(data[, i]), col = 2)
  
  qqnorm(as.numeric(data.rz[, i]),
         main = paste("Normal Q-Q Plot", "RankZ",
                      names(data.rz)[i]))
  qqline(as.numeric(data.rz[, i]), col = 2)
  data.log[, i][!is.finite(data.log[, i])] = NA
  qqnorm(as.numeric(data.log[, i]),
         main = paste("Normal Q-Q Plot", "Log",
                      names(data.log)[i]))
  qqline(data.log[, i], col = 2)
  
}
graphics.off()

#RZ transform seems to work well
