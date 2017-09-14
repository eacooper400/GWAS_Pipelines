### Set the path to wherever you want output to go
### Note that GAPIT will create several output files for each phenotype
setwd("/your-path/")

### Read in the genotype file (in hapmap format)
geno.data=read.table("Input.hmp.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)

### Read in the phenotype file
pheno.data=read.csv("Input.phenotypes.csv", header=TRUE, stringsAsFactors=FALSE)

### GAPIT set up
source("http://www.bioconductor.org/biocLite.R") 
biocLite("multtest")
install.packages("gplots") 
install.packages("LDheatmap") 
install.packages("genetics") 
library("multtest") 
library("gplots") 
library("LDheatmap")
library("genetics") 
library("compiler") #this is in standard R package

source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
source("http://www.zzlab.net/GAPIT/emma.txt")

### Running GAPIT
GAPIT_MLM=GAPIT(Y=pheno.data, G=geno.data)
