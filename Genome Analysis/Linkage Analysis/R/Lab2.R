#read in TDT results from plink
tdt = read.table("C:\\TEMP\\datasets\\plink.tdt", header=T)
dim(tdt)
tdt[1:4,] 
tdt.chr6 <- tdt[tdt$CHR==6,]
dim(tdt.chr6)
#so they are all chromosome 6

# plot
par(mar=c(8,5,5,5))
plot(-log10(tdt.chr6$P), type="n",
     xaxt="n", xlab="", ylab="-log10(p-value)",
     main="Distribution of p-values from TDT across chromosome 6",
     col = "black")
xtick<-seq(1, 1422, by=140)
axis(side=1,at=xtick,labels=tdt.chr6$BP[xtick], las=2)
lines(-log10(tdt.chr6$P),
        type = "h", col = "black")
abline(2.0,0,col="red",lty="dashed")
mtext("Position", side=1, line=6)

plessthan05 <- tdt.chr6[tdt.chr6$P < 0.05,]
dim(plessthan05)

source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")
browseVignettes("qvalue")
library("qvalue")
qobj <- qvalue(p = tdt.chr6$P, fdr.level=0.05)
names(qobj)
qobj$lfdr

p2lessthan05 <- qobj$pvalues[qobj$pvalues < 0.05]
qlessthan05 <- qobj$qvalues[qobj$qvalues < 0.05]
