#read in Fisher's test results from plink output
data = read.table("C:\\TEMP\\datasets\\plink.assoc.fisher", header=T)
dim(data)
data[1:4,]
colnames(data)

# plot
par(mar=c(8,5,5,5))
plot(-log10(data$P), type="n",
     xaxt="n", xlab="", ylab="-log10(p-value)",
     main="Distribution of p-values from Fisher's Test",
     col = "black")
xtick<-seq(1, 1668, by=166)
axis(side=1,at=xtick,labels=data$BP[xtick], las=2)
lines(-log10(data$P),
      type = "h", col = "black")
abline(2.0,0,col="red",lty="dashed")
mtext("Position", side=1, line=6)

plessthan01 <- data[data$P < 0.01,]
dim(plessthan01)
plessthan05 <- data[data$P < 0.05,]
dim(plessthan05)

#read MDS results from plink and plot
mds = read.table("C:\\TEMP\\datasets\\plink.mds", header=T)
colnames(mds)
mds
plot.df <- data.frame(pc1=mds$C1, pc2=mds$C2) 
plot(plot.df, col=c(2,4), xlab="Eigenvector 1",         
     ylab="Eigenvector 2", main="MDS eigenvector 1 vs. eigenvector 2")
legend(0.1, -0.1, c("group 1", "group 2"), col = c(2,4),pch = c(1,1))

mycov <- mds[,c(1,2,4,5)]
write.table(mycov,file="C:\\TEMP\\datasets\\mycov.txt", row.names=FALSE)

covar = read.table("C:\\TEMP\\datasets\\plink.assoc.logistic", header=T)
dim(covar)
colnames(covar)
covar.add <- covar[covar$TEST=="ADD",]
dim(covar.add)

par(mar=c(8,5,5,5))
plot(-log10(covar.add$P), type="n",
     xaxt="n", xlab="", ylab="-log10(p-value)",
     main="Distribution of p-values from Linear Regression",
     col = "black")
xtick<-seq(1, 1668, by=166)
axis(side=1,at=xtick,labels=covar.add$BP[xtick], las=2)
lines(-log10(covar.add$P),
      type = "h", col = "black")
abline(2.0,0,col="red",lty="dashed")
mtext("Position", side=1, line=6)

plessthan01.covar <- covar.add[covar.add$P < 0.01,]
dim(plessthan01.covar)
plessthan05.covar <- covar.add[covar.add$P < 0.05,]
dim(plessthan05.covar)
