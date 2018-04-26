#1
cn = read.table("C:\\TEMP\\datasets\\cn_states.txt", header=T)
logratio = read.table("C:\\TEMP\\datasets\\log2ratios.txt", header=T)
dim(cn)
dim(logratio)

#2 subset the data by chromosome 13
cn.13 <- cn[cn$Chromosome==13,]
logratio.13 <- logratio[logratio$Chromosome==13,]
dim(cn.13)
dim(logratio.13)

# run smoothing and segmentation
source("https://bioconductor.org/biocLite.R")
biocLite("DNAcopy")
library(DNAcopy)
dx <- logratio.13[,-c(1:3)]
d.logratio <- CNA(genomdat=as.matrix(dx),chrom=logratio.13$Chromosome,
                  maploc=logratio.13$Position,data.type=c("logratio"),sampleid=NULL)
d.smoothed <- smooth.CNA(d.logratio)
d.segment <- segment(d.smoothed, verbose=1)

#plot results for chr 1 for all subjects
pdf("CNV_logratio_plot.pdf")
plot(d.segment, plot.type="chrombysample", pt.cex=0.5,lwd=0.5)
dev.off()

#3 plot cn values
cx <- cn.13[,-c(1:3)]
d.cx <- CNA(genomdat=as.matrix(cx),chrom=cn.13$Chromosome,
            maploc=cn.13$Position,data.type=c("logratio"),sampleid=NULL)
d.cx.smoothed <- smooth.CNA(d.cx)
d.cx.segment <- segment(d.cx.smoothed, verbose=1)

# THIS IS WRONG
pdf("CNV_CNStates_plot.pdf")
plot(d.cx.segment, plot.type="chrombysample", main="CN states\nAll Chromosomes")
dev.off()

#CORRECTED CN STATE PLOTS
cx <- cn.13[,-c(1:2)]
cx
pdf("CNV_chr13_cs_state_plot.pdf")
par(mfrow=c(3,3))
for(i in 1:9) {
  plot(range(cx$Position),c(0,4),xlab="Position",ylab="CN State",
       main=paste("Sample",i),type="n")
        points(cx$Position,cx[,(i+1)],col=4,cex=0.35)
        abline(h=2,col=1,lty=2)
}
dev.off()


#4
source("https://bioconductor.org/biocLite.R")
biocLite("CNTools")
library(CNTools)

#create a data matrix of the combined CN vectors
seg <- CNSeg(d.segment$output)
rs.region <- getRS(seg,by="region", imput=FALSE, XY=FALSE,what = "mean")
mat <- rs(rs.region)
mat


#plot the correlation matrix
row.names(mat) <- mat$start
mat.data <- mat[,4:12]
mat.data <- data.matrix(mat.data)
mat.cor <- cor(mat.data,use="pairwise.complete.obs")

par(oma=c(3,3,1,1))
colors <- colorRampPalette(c("red", "white", "blue"))(20)
image(mat.cor,main="Pearson's correlation matrix\nProstate Cancer CNV Samples",axes=F,col=colors)
axis(1,at=seq(0,1,length=ncol(mat.cor)),label=dimnames(mat.cor)[[2]],cex.axis=0.9,las=2)
axis(2,at=seq(0,1,length=ncol(mat.cor)),label=dimnames(mat.cor)[[2]],cex.axis=0.9,las=2)

#5 read in subsetted file
subset = read.table("C:\\TEMP\\datasets\\log2ratio5columns.txt", header=T, row.names=1)
dim(subset)
install.packages("numDeriv")
install.packages(pkgs="c:\\downloads\\ABSOLUTE_1.0.6.tar.gz", repos=NULL, type="source")
library(numDeriv)
library(ABSOLUTE)


