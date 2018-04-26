setwd("~/Datasets/diff")
source("https://bioconductor.org/biocLite.R")
biocLite("cummeRbund")
library(cummeRbund)

#read all cuffdiff files in working directory
cuff <- readCufflinks()
cuff

#density curve
csDensity(genes(cuff))

#find differentially expressed genes
gene_diff_data <- diffData(genes(cuff))
sig_gene_data <- subset(gene_diff_data, (significant == 'yes'))
nrow(sig_gene_data)
sig_gene_data

#print stats for the gene and count table for the isoforms of the gene
tss45_gene <- getGene(cuff, "TSS45")
tss45_gene
head(fpkm(isoforms(tss45_gene)))

# create an expression bar plot for gene isoforms
igb <- expressionBarplot(isoforms(tss45_gene), replicates=T)
igb

# find genes with similar expression patterns
mySimiliar <- findSimilar(cuff,"TSS45", n=4)
mysSimiliar.expression <- expressionPlot(mySimiliar, logMode=T, showErrorbars = F)
mySimiliar
