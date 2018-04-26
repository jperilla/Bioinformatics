#PEARSON's CORRELATION

#create pearson's correlation
install.packages('gplots')
library(gplots)
# method = pearson is the default
data.cdc15cor <- cor(cdc15Data, use="pairwise.complete.obs") # only use pairwise observations where there is data on both side

#plot the correlation matrix
help(heatmap)
rc <- rainbow(nrow(data.cdc15cor), start=0, end=.3)
colors <- colorRampPalette(c("black", "blue", "yellow"))(20)
heatmap(x = data.cdc15cor, col = colors, symm = TRUE, 
        main = "Pearson's correlation - cdc15",
        xlab="cdc15 expression over time", ylab="cdc15 expression over time",
        margins = c(7,7), Rowv=NA, Colv=NA)

install.packages('heatmap3')
library(heatmap3)

heatmap3(x = data.cdc15cor, col = colors, symm = TRUE, 
         main = "Pearson's correlation - cdc15",
         xlab="cdc15 expression over time", ylab="cdc15 expression over time",
         margins = c(8,8), Rowv=NA, Colv=NA)

# PROFILE PLOT
yal002w <- cdc15Data[2,]
yal002w.values <- yal002w[!is.na(yal002w)]
yal002w.mean <- mean(yal002w_values) 
yal002w[is.na(yal002w)] <- yal002w.mean
yal002w

length(yal002w)
plot(range(1:24), range(yal002w), type='n',
     main="YAL002w profile plot", xlab="samples", 
     ylab="intensities")
lines(c(1:24), yal002w)
grid(col='gray')

# BOX PLOT
#pick a gene
gene1.class1 <- data.class1[1,]
gene1.class1.omitna <- gene1.class1[!is.na(gene1.class1)]
gene1.class1.omitna
gene1.class2 <-data.class2[1,]
gene1.class2.omitna <- gene1.class2[!is.na(gene1.class2)]
gene1.class2.omitna

boxplot(gene1.class1.omitna,gene1.class2.omitna, 
        col=c('red','blue'),
        names=c('class 1', 'class 2'),
        main='Example gene from DLBCL cDNA 2-channel dataset',
        axes=F,ylab='log2(ratio intensity)')
axis(2)
axis(1,at=c(1,2),c("class 1","class 2"))

# HISTOGRAM

#plot histograms of both classes
help(par)
par(mfrow=c(2,1)) 
hist(gene1.class1.omitna, 
     main='Example gene from DLBCL cDNA 2-channel dataset - class 1',
     xlab='intensities')

hist(gene1.class2.omitna, 
     main='Example gene from DLBCL cDNA 2-channel dataset - class 2',
     xlab='intensities')


# WELCH's T-TEST
# calculate two-sample Welch’s t-test (unequal variances) between normal and tumor for gene #8000

x <- gene1.class1.omitna
y <- gene1.class2.omitna

# pooled variance
nx <- length(x)
ny <- length(y)
pool.var <- (((nx-1)*var(x)) + ((ny-1)*var(y)))/(nx+ny-2)
pool.var

#minimum sample size calculation
install.packages('pwr')
library(pwr)
dif <- abs(mean(x)-mean(y))/sqrt(pool.var)
dif.1.5fold <- log2(1.5)/sqrt(pool.var)
pl.ss1.5 <- pwr.t.test(d=dif.1.5fold,sig.level=.01,power=0.8,type="two.sample")
pl.ss1.5

pl.ssdelta <- pwr.t.test(d=dif,sig.level=.01,power=0.8,type="two.sample")
pl.ssdelta

#HISTOGRAM OF STANDARD DEVIATIONS

# multiple gene power curves using the help file data set
install.packages('samplesize')
library(samplesize)
library(gdata)
data.withsd <- transform(data, SD=apply(data,1, sd, na.rm = TRUE))
dim(data.withsd)
data.sd <- data.withsd[,41]
hist(data.sd, main='Histogram of standard deviations for 13,412 genes',
     ylab='Frequency', xlab='Standard deviation (for data on the log scale)')


# SAMPLE SIZE PLOT

biocLite("ssize")
library(ssize)

# plot a gene proportion vs. sample size plot using the same criteria AND power=80%
all.size <- ssize(sd=data.sd, delta=log2(3.0), sig.level=0.05, power=0.8) 
ssize.plot(all.size, lwd=2, col="magenta", xlim=c(1,20)) 
xmax <- par("usr")[2]-1; 
ymin <- par("usr")[3] + 0.05 
legend(x=xmax, y=ymin, legend= strsplit( paste("fold change=",3.0,",", "alpha=", 0.05, ",", "power=",0.8,",", "# genes=", length(data.sd), sep=''), "," )[[1]], xjust=1, yjust=0, cex=1.0) 
title("Sample Size to Detect 2-Fold Change") 
