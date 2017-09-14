# 17 Dec. 2013 A.A. Comeault
# Script for importing, calculating the mean, and analyzing output from gemma BSLMM. 
# This file combines runs of probit BSLMM carried out WITH 10 CHAINS.

# select the appropriate wd containing the gemma output files.
setwd("path_to_output_files/")


rm(list=ls())
library(ggplot2)

##########################
# hyperparameter summary #
##########################
# 1) import and calculate the mean .hyp.txt from 10 chains of BSLMM run in gemma
# this can take a long time. on my computer it took about 10 or 15 minutes.
hyp.files <- list.files(pattern = "*.hyp.txt")

for (i in seq_along(hyp.files)) {
  if (i == 1) {
    hyp.dat <- read.table(hyp.files[i], header=TRUE)
  } else {
    dat <- read.table(hyp.files[i], header=TRUE)
    hyp.dat <- cbind(hyp.dat,dat)
  }
}
# 2) calculate the mean hyperparameter estimates across the 10 chains
# again, takes a while.
hyp.mean <- data.frame(cbind(apply(hyp.dat[,seq(1,dim(hyp.dat)[2],6)], MARGIN=1, FUN=mean, na.rm=TRUE),
                             apply(hyp.dat[,seq(2,dim(hyp.dat)[2],6)], MARGIN=1, FUN=mean, na.rm=TRUE),
                             apply(hyp.dat[,seq(3,dim(hyp.dat)[2],6)], MARGIN=1, FUN=mean, na.rm=TRUE),
                             apply(hyp.dat[,seq(4,dim(hyp.dat)[2],6)], MARGIN=1, FUN=mean, na.rm=TRUE),
                             apply(hyp.dat[,seq(5,dim(hyp.dat)[2],6)], MARGIN=1, FUN=mean, na.rm=TRUE),
                             apply(hyp.dat[,seq(6,dim(hyp.dat)[2],6)], MARGIN=1, FUN=mean, na.rm=TRUE)))
colnames(hyp.mean) <- c("h","pve","rho","pge","pi","n_gamma")
rm(list=ls(pattern="dat"))

# 3) Generate the table for median and 95% EPTI estimates for h, PVE, rho, PGE
# 50% quantile == median values, 0.025 and 0.975 are lower and upper 95%ETPI
h     <- quantile(hyp.mean$h, probs = c(0.5, 0.025, 0.975))
PVE   <- quantile(hyp.mean$pve, probs = c(0.5, 0.025, 0.975))
rho   <- quantile(hyp.mean$rho, probs = c(0.5, 0.025, 0.975))
PGE   <- quantile(hyp.mean$pge, probs = c(0.5, 0.025, 0.975))
n_gam <- quantile(hyp.mean$n_gamma, probs = c(0.5,0.025, 0.975))

table <- rbind(h,rho,PVE,PGE,n_gam)
colnames(table) <- c("median", "low95EPTI", "upp95EPTI")

# save table:
write.table(table,"Table_file_name.csv",quote=FALSE,sep=",",row.names=TRUE,col.names=TRUE)


##########################
# individual SNP effects #
##########################

# 1) import and calculate the mean .param.txt from 10 chains of BSLMM run in gemma 
para.files <- list.files(pattern = "*.param.txt")

for (i in seq_along(para.files)) {
  if (i == 1) {
    para.dat <- read.table(para.files[i], header=TRUE)
  } else {
    dat <- read.table(para.files[i], header=TRUE)
    para.dat <- cbind(para.dat, dat[,4:7])
  }
}

# 2) calculate mean alpha, beta, and gamma across all 10 chains
para.mean <- data.frame(cbind(para.dat[,1:3],
                              apply(para.dat[,seq(5,dim(para.dat)[2],4)], MARGIN=1, FUN=mean),
                              apply(para.dat[,seq(6,dim(para.dat)[2],4)], MARGIN=1, FUN=mean),
                              apply(para.dat[,seq(7,dim(para.dat)[2],4)], MARGIN=1, FUN=mean)))
colnames(para.mean) <- c("CHR","RS","PS","alpha","beta","gamma")
para.mean$abs.beta.g  <- abs(para.mean$beta * para.mean$gamma) # calculate beta x gamma
rm(list=ls(pattern="dat"))

# 3) make subset datasets containing snps satisfying certain conditions:
sparse.effects <- para.mean[para.mean$abs.beta.g != 0,] # remove snps with no sparse effect on phenotypes
top1snp.effects   <- sparse.effects[which(sparse.effects$abs.beta.g > quantile(x=sparse.effects$abs.beta.g,probs=0.99)),] # get top 1% high effect snps
top0.1snp.effects  <- sparse.effects[which(sparse.effects$abs.beta.g > quantile(x=sparse.effects$abs.beta.g,probs=0.999)),] # get top 0.1% high effect snps
pip01       <- para.mean[which(para.mean$gamma >= 0.01),] # snps with gamma (i.e. PIP) > 0.01
pip10       <- para.mean[which(para.mean$gamma >= 0.40),] # gamma > 0.40
pip50       <- para.mean[which(para.mean$gamma >= 0.25),] # gamma > 0.25

# Save what you'd like to:
write.table(pip01, "File_name_here.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(para.mean, "File_name_here.txt", quote=FALSE, row.names=FALSE, sep="\t")
