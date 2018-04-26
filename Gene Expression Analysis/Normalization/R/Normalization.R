source("https://bioconductor.org/biocLite.R")
biocLite("marray")
a
library(marray)

# extract raw marray data from genepix
help(read.GenePix)
GSM304445 <- read.GenePix("GSM304445.gpr", "c:\\temp\\datasets\\")
GSM304446 <- read.GenePix("GSM304446.gpr", "c:\\temp\\datasets\\")
GSM304447 <- read.GenePix("GSM304447.gpr", "c:\\temp\\datasets\\")
GSM304448 <- read.GenePix("GSM304448.gpr", "c:\\temp\\datasets\\")

class?marrayRaw
GSM304445.Cy5foreground.median <- median(targets.GSM304445@maRf)
GSM304445.Cy5background.median <- median(targets.GSM304445@maRb)
GSM304445.Cy3foreground.median <- median(targets.GSM304445@maGf)
GSM304445.Cy3background.median <- median(targets.GSM304445@maGb)

GSM304446.Cy5foreground.median <- median(targets.GSM304446@maRf)
GSM304446.Cy5background.median <- median(targets.GSM304446@maRb)
GSM304446.Cy3foreground.median <- median(targets.GSM304446@maGf)
GSM304446.Cy3background.median <- median(targets.GSM304446@maGb)

GSM304447.Cy5foreground.median <- median(targets.GSM304447@maRf)
GSM304447.Cy5background.median <- median(targets.GSM304447@maRb)
GSM304447.Cy3foreground.median <- median(targets.GSM304447@maGf)
GSM304447.Cy3background.median <- median(targets.GSM304447@maGb)

GSM304448.Cy5foreground.median <- median(targets.GSM304448@maRf)
GSM304448.Cy5background.median <- median(targets.GSM304448@maRb)
GSM304448.Cy3foreground.median <- median(targets.GSM304448@maGf)
GSM304448.Cy3background.median <- median(targets.GSM304448@maGb)

#3 types of normalization
GSM304445 <- targets.GSM304445
GSM304445.norm.median <- maNorm(GSM304445, norm="median")
GSM304445.norm.loess <- maNorm(GSM304445, norm="loess")
GSM304445.norm.print <- maNorm(GSM304445, norm="scalePrintTipMAD")

#plots
par(mfrow=c(2,2))
plot(GSM304445, lines.func=NULL, legend=NULL,
     main="GSM304445: pre-normalization MA--plot")
plot(GSM304445.norm.median, lines.func=NULL, legend=NULL,
     main="GSM304445: median global normalization MA--plot")
plot(GSM304445.norm.loess, lines.func=NULL, legend=NULL,
     main="GSM304445: loess normalization MA--plot")
plot(GSM304445.norm.print, lines.func=NULL, legend=NULL,
     main="GSM304445: print tip normalization MA--plot")

#GSM304446 -3 types of normalization
GSM304446 <- targets.GSM304446
GSM304446.norm.median <- maNorm(GSM304446, norm="median")
GSM304446.norm.loess <- maNorm(GSM304446, norm="loess")
GSM304446.norm.print <- maNorm(GSM304446, norm="scalePrintTipMAD")

#GSM304446 - plots
par(mfrow=c(2,2))
plot(GSM304446, lines.func=NULL, legend=NULL,
     main="GSM304446: pre-normalization MA--plot")
plot(GSM304446.norm.median, lines.func=NULL, legend=NULL,
     main="GSM304446: median global normalization MA--plot")
plot(GSM304446.norm.loess, lines.func=NULL, legend=NULL,
     main="GSM304446: loess normalization MA--plot")
plot(GSM304446.norm.print, lines.func=NULL, legend=NULL,
     main="GSM304446: print tip normalization MA--plot")

#GSM304447 -3 types of normalization
GSM304447 <- targets.GSM304447
GSM304447.norm.median <- maNorm(GSM304447, norm="median")
GSM304447.norm.loess <- maNorm(GSM304447, norm="loess")
GSM304447.norm.print <- maNorm(GSM304447, norm="scalePrintTipMAD")

#GSM304447 - plots
par(mfrow=c(2,2))
plot(GSM304447, lines.func=NULL, legend=NULL,
     main="GSM304447: pre-normalization MA--plot")
plot(GSM304447.norm.median, lines.func=NULL, legend=NULL,
     main="GSM304447: median global normalization MA--plot")
plot(GSM304447.norm.loess, lines.func=NULL, legend=NULL,
     main="GSM304447: loess normalization MA--plot")
plot(GSM304447.norm.print, lines.func=NULL, legend=NULL,
     main="GSM304447: print tip normalization MA--plot")

#GSM304448 -3 types of normalization
GSM304448 <- targets.GSM304448
GSM304448.norm.median <- maNorm(GSM304448, norm="median")
GSM304448.norm.loess <- maNorm(GSM304448, norm="loess")
GSM304448.norm.print <- maNorm(GSM304448, norm="scalePrintTipMAD")

#GSM304448 - plots
par(mfrow=c(2,2))
plot(GSM304448, lines.func=NULL, legend=NULL,
     main="GSM304448: pre-normalization MA--plot")
plot(GSM304448.norm.median, lines.func=NULL, legend=NULL,
     main="GSM304448: median global normalization MA--plot")
plot(GSM304448.norm.loess, lines.func=NULL, legend=NULL,
     main="GSM304448: loess normalization MA--plot")
plot(GSM304448.norm.print, lines.func=NULL, legend=NULL,
     main="GSM304448: print tip normalization MA--plot")


#density plot of log rations
par(mfrow=c(1,1))
density <- density(maM(GSM304448.norm.loess), na.rm=TRUE)
plot(x=NULL, y=NULL,
     xlim=c(min(density$x), max(density$x)),
     ylim=c(min(density$y), max(density$y)), type='n',
     main="Density plots of log-ratios M",
     xlab="Log ratios", ylab="Density")
lines(density(maM(GSM304448), na.rm=TRUE), lwd=2, col=1)
lines(density(maM(GSM304448.norm.median), na.rm=TRUE), lwd=2, col=4 )
lines(density(maM(GSM304448.norm.loess), na.rm=TRUE), lwd=2, col=3 )
lines(density(maM(GSM304448.norm.print), na.rm=TRUE), lwd=2, col=2 )
help(legend)
legend(x=2.5, y=0.7, legend=c("None", "Median", "Lowess", "Print-tip Lowess"),
       col=c(1, 4, 3, 2), pch = c(NA, 3, 4), lty=1, lwd=2)

#single-channel normalization
help(limma)
limmaUsersGuide(view=TRUE)

GSM304445.Cy5.f <- GSM304445@maRf
GSM304445.Cy5.b <- GSM304445@maRb
GSM304445.Cy5.subtract <- GSM304445.Cy5.f - GSM304445.Cy5.b
GSM304445.Cy5.subtract.log <- log2(GSM304445.Cy5.subtract)
GSM304445.Cy5.subtract.log
GSM304445.Cy5.norm <- normalizeBetweenArrays(GSM304445.Cy5.subtract.log, method="scale")

GSM304446.Cy5.f <- GSM304446@maRf
GSM304446.Cy5.b <- GSM304446@maRb
GSM304446.Cy5.subtract <- GSM304446.Cy5.f - GSM304446.Cy5.b
GSM304446.Cy5.subtract.log <- log2(GSM304446.Cy5.subtract)
GSM304446.Cy5.subtract.log
GSM304446.Cy5.norm <- normalizeBetweenArrays(GSM304446.Cy5.subtract.log, method="scale")

GSM304447.Cy5.f <- GSM304447@maRf
GSM304447.Cy5.b <- GSM304447@maRb
GSM304447.Cy5.subtract <- GSM304447.Cy5.f - GSM304447.Cy5.b
GSM304447.Cy5.subtract.log <- log2(GSM304447.Cy5.subtract)
GSM304447.Cy5.subtract.log
GSM304447.Cy5.norm <- normalizeBetweenArrays(GSM304447.Cy5.subtract.log, method="scale")

GSM304448.Cy5.f <- GSM304448@maRf
GSM304448.Cy5.b <- GSM304448@maRb
GSM304448.Cy5.subtract <- GSM304448.Cy5.f - GSM304448.Cy5.b
GSM304448.Cy5.subtract.log <- log2(GSM304448.Cy5.subtract)
GSM304448.Cy5.subtract.log
GSM304448.Cy5.norm <- normalizeBetweenArrays(GSM304448.Cy5.subtract.log, method="scale")

help(cor.test)
GSM56.cor <- cor(as.numeric(GSM304445.Cy5.norm), as.numeric(GSM304446.Cy5.norm), 
                 method = "spearman", use="pairwise.complete.obs")
GSM57.cor <- cor(as.numeric(GSM304445.Cy5.norm), as.numeric(GSM304447.Cy5.norm), 
                 method = "spearman", use="pairwise.complete.obs")
GSM58.cor <- cor(as.numeric(GSM304445.Cy5.norm), as.numeric(GSM304448.Cy5.norm), 
                 method = "spearman", use="pairwise.complete.obs")
GSM67.cor <- cor(as.numeric(GSM304446.Cy5.norm), as.numeric(GSM304447.Cy5.norm), 
                 method = "spearman", use="pairwise.complete.obs")
GSM68.cor <- cor(as.numeric(GSM304446.Cy5.norm), as.numeric(GSM304448.Cy5.norm), 
                 method = "spearman", use="pairwise.complete.obs")
GSM78.cor <- cor(as.numeric(GSM304447.Cy5.norm), as.numeric(GSM304448.Cy5.norm), 
                 method = "spearman", use="pairwise.complete.obs")

GSM56.loess.cor <- cor(maM(GSM304445.norm.loess), maM(GSM304446.norm.loess), 
                 method = "spearman", use="pairwise.complete.obs")
GSM57.loess.cor <- cor(maM(GSM304445.norm.loess), maM(GSM304447.norm.loess), 
                       method = "spearman", use="pairwise.complete.obs")
GSM58.loess.cor <- cor(maM(GSM304445.norm.loess), maM(GSM304448.norm.loess), 
                       method = "spearman", use="pairwise.complete.obs")
GSM67.loess.cor <- cor(maM(GSM304446.norm.loess), maM(GSM304447.norm.loess), 
                       method = "spearman", use="pairwise.complete.obs")
GSM68.loess.cor <- cor(maM(GSM304446.norm.loess), maM(GSM304448.norm.loess), 
                       method = "spearman", use="pairwise.complete.obs")
GSM78.loess.cor <- cor(maM(GSM304447.norm.loess), maM(GSM304448.norm.loess), 
                       method = "spearman", use="pairwise.complete.obs")

help(pairs)
Cy5.norm.matrix <- cbind(GSM304445.Cy5.norm, GSM304446.Cy5.norm,
                         GSM304447.Cy5.norm, GSM304448.Cy5.norm)
Cy5.norm.matrix[1:4, 1:4]
colnames(Cy5.norm.matrix) <- c("GSM30445", "GSM30446", "GSM30447", "GSM30448")
loess.norm.matrix <- cbind(maM(GSM304445.norm.loess), maM(GSM304446.norm.loess),
                         maM(GSM304447.norm.loess), maM(GSM304448.norm.loess))
colnames(loess.norm.matrix) <- c("GSM30445", "GSM30446", "GSM30447", "GSM30448")

loess.norm.matrix[1:4, 1:4]
pairs(Cy5.norm.matrix, labels=colnames(Cy5.norm.matrix),
      main="Cy5 single channel Scatterplot Matrix")

pairs(Cy5.norm.matrix, labels=colnames(Cy5.norm.matrix),
      main="Loess Normalization Scatterplot matrix")

cor.norm.table <- rbind(c(1.0, GSM56.loess.cor, GSM57.loess.cor, GSM58.loess.cor),
                   c(GSM56.loess.cor, 1.0, GSM67.loess.cor, GSM68.loess.cor),
                   c(GSM57.loess.cor, GSM67.loess.cor, 1.0, GSM78.loess.cor),
                   c(GSM58.loess.cor, GSM68.loess.cor, GSM78.loess.cor, 1.0))
rownames(cor.norm.table) <- c("GSM304445", "GSM304446", "GSM304447", "GSM304448")
colnames(cor.norm.table) <- c("GSM304445", "GSM304446", "GSM304447", "GSM304448")
cor.norm.table

cor.norm.table <- rbind(c(1.0, GSM56.cor, GSM57.cor, GSM58.cor),
                       c(GSM56.cor, 1.0, GSM67.cor, GSM68.cor),
                       c(GSM57.cor, GSM67.cor, 1.0, GSM78.cor),
                       c(GSM58.cor, GSM68.cor, GSM78.cor, 1.0))
rownames(cor.Cy5.table) <- c("GSM304445", "GSM304446", "GSM304447", "GSM304448")
colnames(cor.Cy5.table) <- c("GSM304445", "GSM304446", "GSM304447", "GSM304448")
cor.Cy5.table




GSM304445.sub <- GSM304445@maRf - GSM304445@maRb
GSM304446.sub <- GSM304446@maRf - GSM304446@maRb
GSM304447.sub <- GSM304447@maRf - GSM304447@maRb
GSM304448.sub <- GSM304448@maRf - GSM304448@maRb
GSM304445.sorted <- sort(GSM304445.sub)
GSM304446.sorted <- sort(GSM304446.sub)
GSM304447.sorted <- sort(GSM304447.sub)
GSM304448.sorted <- sort(GSM304448.sub)

GSM.matrix <- cbind(GSM304445.sorted, GSM304446.sorted, 
                    GSM304447.sorted, GSM304448.sorted)

rowmeans <- rowMeans(GSM.matrix, na.rm = FALSE, dims = 1)
GSM.matrix <- cbind(GSM.matrix, rowmeans)
GSM304445.rank <- rank(GSM304445.sub, ties="first")
GSM304446.rank <- rank(GSM304446.sub, ties="first")
GSM304447.rank <- rank(GSM304447.sub, ties="first")
GSM304448.rank <- rank(GSM304448.sub, ties="first")

help(match)
GSM.matrix <- cbind(GSM.matrix, GSM304445.rank, 
                    GSM304446.rank, GSM304447.rank, 
                    GSM304448.rank)

help(sort)
GSM.matrix[1:4,]
GSM.matrix[,1]<- GSM.matrix[order(GSM304445.rank),][,1]
GSM.matrix[,2]<- GSM.matrix[order(GSM304446.rank),][,2]
GSM.matrix[,3]<- GSM.matrix[order(GSM304447.rank),][,3]
GSM.matrix[,4]<- GSM.matrix[order(GSM304448.rank),][,4]
GSM.matrix[1:4,]

par(mfrow=c(2,2))
hist(GSM.matrix[,1])
hist(GSM.matrix[,2])
hist(GSM.matrix[,3])
hist(GSM.matrix[,4])

GSM.matrix.log <- log2(GSM.matrix[,1:4])
GSM.matrix.log[1:4,]

GSM56.cor <- cor(as.numeric(GSM.matrix.log[,1]), as.numeric(GSM.matrix.log[,2]), 
                 method = "spearman", use="pairwise.complete.obs")
GSM57.cor <- cor(as.numeric(GSM.matrix.log[,1]), as.numeric(GSM.matrix.log[,3]),  
                 method = "spearman", use="pairwise.complete.obs")
GSM58.cor <- cor(as.numeric(GSM.matrix.log[,1]), as.numeric(GSM.matrix.log[,4]),  
                 method = "spearman", use="pairwise.complete.obs")
GSM67.cor <- cor(as.numeric(GSM.matrix.log[,2]), as.numeric(GSM.matrix.log[,3]),  
                 method = "spearman", use="pairwise.complete.obs")
GSM68.cor <- cor(as.numeric(GSM.matrix.log[,2]), as.numeric(GSM.matrix.log[,4]), 
                 method = "spearman", use="pairwise.complete.obs")
GSM78.cor <- cor(as.numeric(GSM.matrix.log[,3]), as.numeric(GSM.matrix.log[,4]),  
                 method = "spearman", use="pairwise.complete.obs")
""
colnames(GSM.matrix.log) <- c("GSM304445", "GSM304446", "GSM304447", "GSM304448")
pairs(GSM.matrix.log, labels=colnames(GSM.matrix.log),
      main="Quantile Normalized Scatterplot Matrix")

quantile.norm.table <- rbind(c(1.0, GSM56.cor, GSM57.cor, GSM58.cor),
                        c(GSM56.cor, 1.0, GSM67.cor, GSM68.cor),
                        c(GSM57.cor, GSM67.cor, 1.0, GSM78.cor),
                        c(GSM58.cor, GSM68.cor, GSM78.cor, 1.0))
quantile.norm.table


# qt-PCR
f.parse <- function(path=pa,file=fi,out=out.fi) {
  d <- read.table(paste(path,file,sep=""),skip=11,sep=",",header=T)
  u <- as.character(unique(d$Name))
  u <- u[u!=""]; u <- u[!is.na(u)];
  ref <- unique(as.character(d$Name[d$Type=="Reference"]))
  u <- unique(c(ref,u))
  hg <- c("B-ACTIN","GAPDH","18S")
  hg <- toupper(hg)
  p <- unique(toupper(as.character(d$Name.1)))
  p <- sort(setdiff(p,c("",hg)))
  
  mat <- matrix(0,nrow=length(u),ncol=length(p))
  dimnames(mat) <- list(u,p)
  for (i in 1:length(u)) {
    print(paste(i,": ",u[i],sep=""))
    tmp <- d[d$Name %in% u[i],c(1:3,6,9)]
    g <- toupper(unique(as.character(tmp$Name.1)))
    g <- sort(setdiff(g,c("",hg)))
    
    for (j in 1:length(g)) {
      v <- tmp[toupper(as.character(tmp$Name.1)) %in% g[j],5]
      v <- v[v!=999]
      v <- v[((v/mean(v))<1.5) & ((v/mean(v))>0.67)]	#gene j vector
      
      hv3 <- NULL
      for (k in 1:length(hg)) {	#housekeeping gene vector (each filtered by reps)
        hv <- tmp[toupper(as.character(tmp$Name.1)) %in% hg[k],5]
        hv <- hv[hv!=999]
        hv3 <- c(hv3,hv[((hv/mean(hv))<1.5) & ((hv/mean(hv))>0.67)]) 	
      }
      sv <- mean(as.numeric(v)) - mean(as.numeric(hv3))	#scaled value for gene j
      
      if(i==1) { #reference sample only
        mat[u[i],g[j]] <- sv
        next
      }
      
      mat[u[i],g[j]] <- sv - mat[u[1],g[j]]
    }
  }
  
  mat[1,][!is.na(mat[1,])] <- 0
  fc <- 2^(-1 * mat)
  write.table(t(c("Subject",dimnames(mat)[[2]])),paste(path,out,sep=""),quote=F,sep="\t",col.names=F,row.names=F)
  write.table(round(fc,3),paste(path,out,sep=""),quote=F,sep="\t",append=T,col.names=F)
}

# run function
pa <- "C:\\temp\\datasets\\"
fi <- "qRT-PCR.CSV"
output.fi <- "foldchange.txt"

f.parse(pa,fi,output.fi)


data = read.table("c:\\temp\\datasets\\foldchange.txt", header=T, row.names=1)

