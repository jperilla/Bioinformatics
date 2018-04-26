#get thyroid and testes expression datasets
data.thyroid = read.table("//Volumes/Drive/Bioinformatics/s1.genes.results", header=T)
dim(data.thyroid)
colnames(data.thyroid)
data.thyroid[1:4,]

data.testes = read.table("//Volumes/Drive/Bioinformatics/s2.genes.results", header=T)
dim(data.testes)
colnames(data.testes)
data.testes[1:4,]

#subset for just TPM and FPKM
data.thyroid.stats <- data.thyroid[,6:7]
data.thyroid.stats[1:4,]
data.testes.stats <- data.testes[,6:7]
data.testes.stats[1:4,]

#normalize TPM and FPKM
normalize <- function(x) {
  return(log2(x+1))
}
data.thyroid.stats$TPM <- lapply(data.thyroid.stats$TPM, function(x) sapply(x, normalize))
data.thyroid.stats$FPKM <- lapply(data.thyroid.stats$FPKM, function(x) sapply(x, normalize))
data.testes.stats$TPM <- lapply(data.testes.stats$TPM, function(x) sapply(x, normalize))
data.testes.stats$FPKM <- lapply(data.testes.stats$FPKM, function(x) sapply(x, normalize))

# compare
data.thyroid[1:4,6:7]
data.thyroid.stats[1:4,]
data.testes[1:4,6:7]
data.thyroid.stats[1:4,]

#Pearson's correlation and scatter plot
help("cor")
data.FPKM.cor <- cor(as.numeric(data.thyroid.stats$FPKM), y=as.numeric(data.testes.stats$FPKM), use="pairwise.complete.obs")
data.TPM.cor <- cor(as.numeric(data.thyroid.stats$TPM), y=as.numeric(data.testes.stats$TPM), use="pairwise.complete.obs")
data.FPKM.cor
data.TPM.cor

plot(data.thyroid.stats$FPKM, data.testes.stats$FPKM,
     xlab="Thyroid FPKM values", ylab="Testes FPKM values",
     main="Thyroid vs. Testes FPKM values\nPearson's coefficient=0.8254891",
     col="blue")


plot(data.thyroid.stats$TPM, data.testes.stats$TPM,
     xlab="Thyroid TPM values", ylab="Testes TPM values",
     main="Thyroid vs. Testes TPM values\nPearson's coefficient=0.8256516",
     col="darkgreen")

# calculate fold change between testes and thyroid
foldchange <- as.numeric(data.thyroid.stats$TPM) - as.numeric(data.testes.stats$TPM)
data <- cbind(data.thyroid[,1:2], foldchange)
data.sorted <- data[order(foldchange),] 
dim(data.sorted)
genes.testes.over <- data.sorted[1:10,]$gene_id
genes.thyroid.over <- data.sorted[26100:26109,]$transcript_id
genes.testes.over
genes.thyroid.over

