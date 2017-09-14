### Set the path to the directory with your output
setwd("/your-path")

### Reading in GAPIT Files
gapit.files=list.files(path="your-path-to-GapitFiles", pattern="*GWAS.Results.csv", full.names=TRUE)

### Read in and get the top 1% of the GAPIT results
gapit.top1=list()
for (i in 1:(length(gapit.files)) {
    results=read.csv(gapit.files[i], header=TRUE, stringsAsFactors=FALSE)
    results$logp=(-1) * log10(results$FDR_Adjusted_P.values)
    top1snp=results[which(results$logp > quantile(x=results$logp,probs=0.99)), c(1,2,3,10)] # This gets the top 1%.  If you want something different, change the value of probs=...
    gapit.top1[[i]]=top1snp
}
names(gapit.top1) = c("Names of your traits")

### Do the same for the GEMMA results
gemma.top1=list()
gemma.files=list.files(path="GEMMA_Results", pattern="*top1.snps.txt", full.names=TRUE) # Notice here I had already output only the top 1% of SNPs, so I don't filter anymore...
for (i in 1:(length(gemma.files)) {
    results=read.table(gemma.files[i], header=TRUE, stringsAsFactors=FALSE)
    results=results[order(results$gamma, decreasing=TRUE),]
    gemma.top1[[i]]=results[,c(1,2,3,6)]
}
names(gemma.top1)=c("Names of your traits")

### Read in the gff3 annotation file
gene.file=read.table("your-genome-annotation.gff3", header=FALSE, stringsAsFactors=FALSE)
genes.only=subset(gene.file, gene.file$V3=="gene")

### Function to get all genes in a specified region around a SNP
### By default, this uses a window size of 50kb; you can change the value of regionSize for something different if you want
get.genes=function(chr, pos, regionSize=50000, gene.table=genes.only) {
    g1=gene.table[grep(chr, gene.table$V1),]
    g=subset(g1, g1$V4>=(pos-regionSize) & g1$V5 <= (pos+regionSize))
    name=unlist(lapply(strsplit(g$V9, split=";"), `[[`, 2))
    name=gsub("Name=", "", name)
    return(name)
}

### Here is an example showing how I could use the function
### with the top 1% snp file for a certain trait
my.best.snps = gemma.top1[["My Trait Name"]]
my.best.snps$Genes_in_50kb = rep(NA, nrow(my.best.snps)) # This is just setting up an empty column for results

### Run the function with a 50kb window on either side
for (i in 1:nrow(my.best.snps)) {
    my.genes=paste(get.genes(chr=my.best.snps[i,1], pos=my.best.snps[i,3], regionSize=50000, gene.table=genes.only), collapse=",") # Note the columns for CHR and PS in GEMMA output might be different in GAPIT!!
    my.best.snps$Genes_in_50kb[i]=my.genes
}

### Optional: Get info from the file with gene descriptions
gene.desc=read.table("your path to annotation_info.txt", header=FALSE, stringsAsFactors=FALSE, fill=TRUE)
cand.genes=unlist(strsplit(paste(my.best.snps$Genes_in_50kb, collapse=","), split=","))
m=match(cand.genes, gene.desc$V2)

relevant.gene.desc=gene.desc[m,]

write.table(my.best.snps, "Top_SNPs_GEMMA_w_ClosesGeneIDs.txt", quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(relevant.gene.desc, "Genes_within50kb_topSNPs.txt", quote=FALSE, col.names=TRUE, row.names=FALSE)
