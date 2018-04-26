#Load libraries
source("http://www.bioconductor.org/biocLite.R")
biocLite("GEOquery")
library(Biobase)
library(GEOquery)
help(getGEO)

#Read in dataset, analyze format of data, and set row names to Gene Ids
metadata <- getGEO(filename='C:\\TEMP\\datasets\\final\\GDS4222.soft')
metadata@header
metadata@dataTable # CLASSES ARE IN HERE
metadata@gpl
data <- Table(metadata)
dim(data)
colnames(data)
rownames(data)
data[1:4,1:4]
rownames(data) <- data[,1]
rownames(data)
data.GeneIds <- data[,2]
data <- data[,3:132]

#normalize data, save linear data first
data.linear <- data
data <- log2(data)


# Identify potential outliers

# Correlation Plot Heatmap
install.packages('heatmap3') 
library(heatmap3) 
heatmap3(x = data, symm = FALSE,           
         main = "Correlation Plot",        
         ylab="54,675 probesets", labRow = "",         
         Rowv=NA, Colv=NA, useRaster=TRUE,  cexCol=0.5)

#Hierarchical Clustering Dendogram
data.transpose <- t(data) 
data.transpose.dist <- dist(data.transpose,method="euclidean") 
data.transpose.clust <- hclust(data.transpose.dist,method="single") 
plot(data.transpose.clust,labels=names(data.transpose),cex=0.75, 
     main="Outliers - Cluster Dendogram") 

#CV vs. Mean Plot
data.mean <- apply(data,2,mean) 
data.sd <- sqrt(apply(data,2,var)) 
data.cv <- data.sd/data.mean 
plot(data.mean,data.cv,      
     main="All Samples\n CV vs. Mean",      
     xlab="Mean",ylab="CV",col='blue',      
     cex=1.5,type="n")
points(data.mean,data.cv,bg="lightblue",col=1,pch=21) 
text(data.mean,data.cv,label=dimnames(data)[[2]],pos=3,cex=.8)

#average correlation plot
install.packages('gplots') 
library(gplots) 
data.cor <- cor(data) 
data.avg <- apply(data.cor,1,mean) 
par(oma=c(3,0.1,0.1,0.1)) 
plot(c(1,length(data.avg)),range(data.avg),      
     type="n",xlab="",ylab="Avg r",      
     main="Avg correlation of All samples",      
     axes=F) 
points(data.avg,bg="red",col=1,pch=21,cex=1.25) 
axis(1,at=c(1:length(data.avg)),labels=dimnames(data)[[2]],las=2, cex.lab=0.4,cex.axis=0.6) 
axis(2) 
abline(v=seq(0.5,130.5,1),col="grey")

#Remove outliers
dim(data)
remove1 <- grep("GSM447646", names(data))
remove2 <- grep("GSM447628", names(data))
remove3 <- grep("GSM447655", names(data))
cols.removed <- c(remove1, remove2, remove3)
data <- data[ , -which(names(data) %in% c("GSM447646","GSM447628","GSM447655"))]
data.linear <- data.linear[ , -which(names(data.linear) %in% c("GSM447646","GSM447628","GSM447655"))]
dim(data)
dim(data.linear)

#remove outliers from class too
class <- data.frame(metadata@dataTable@columns[5]) 
dim(class)
class <- class[-cols.removed,,drop=FALSE]
rownames(class) <- NULL
class

#Filter out lowest 5% of expression values
data.genes.mean <- apply(data.linear,1,mean) 
data.genes.mean.sorted <- sort(data.genes.mean)
fivepercent <- round(0.05*length(data.genes.mean))
genestofilter <- names(data.genes.mean.sorted)[1:fivepercent]
data <- data[-which(rownames(data) %in% genestofilter),]
data.linear <- data.linear[-which(rownames(data.linear) %in% genestofilter),]
dim(data)
dim(data.linear)

#Feature Selection
#get the classes, two factor:
# relapse or not (combine all types of relapse into one group: EARLY, LATE, REFRACTORY)

na <- grep("relapse type: NA", class[,1])
e <- grep("relapse type: EARLY", class[,1])
l <- grep("relapse type: LATE", class[,1])
r <- grep("relapse type: REFRACTORY", class[,1])
gc <- na# population that was cured by first line treatment
act <- c(e,l,r)# population that relapsed or had refractory disease
length(gc)
length(act)

# function to calculate Student’s two-sample t-test on all genes at once
# function returns the p-value for the test
# NAs are removed for each test
t.test.all.genes <- function(x,s1,s2) {
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- t.test(x1,x2, alternative="two.sided",var.equal=T)
  out <- as.numeric(t.out$p.value)
  return(out)
}

# distribution of p-values and fold change
pv <- apply(data,1,t.test.all.genes,s1=gc,s2=act)
par(mfrow=c(1,2))
hist(pv,col="lightblue",xlab="p-values",main="P-value dist’n between\nGC and ACT groups",cex.main=0.9)
threshold1 <- 0.05
threshold2 <- 0.01
threshold3 <- 0.0001 #needed to go this low to get a reasonable # of genes
abline(v=threshold1,col=2,lwd=2)
abline(v=threshold2,col=4, lwd=2)
abline(v=threshold3,col=3, lwd=2)
hist(-log10(pv),col="lightblue",xlab="log10(p-values)", main="-log10(pv) dist’n between\nGC and ACT groups",cex.main=0.9)
abline(v= -log10(threshold1),col=2,lwd=2)
abline(v= -log10(threshold2),col=4,lwd=2)
abline(v= -log10(threshold2),col=3,lwd=2)

pv.lessthan05 <- pv[pv < threshold1]
length(pv.lessthan05)
pv.lessthan01 <- pv[pv < threshold2]
length(pv.lessthan01)
pv.sig <- pv[pv < threshold3]
length(pv.sig)
names(pv.sig)
hist(pv.sig, col="lightblue", xlab="p-values", main="Histogram of genes with significant p-values")

#fold change success = no relapse, failure = relapse or refractory
success.m <- apply(data[,gc],1,mean,na.rm=T) 
failure.m <- apply(data[,act],1,mean,na.rm=T) 
fold <- success.m-failure.m
summary(fold)
fold.sig <- fold[abs(fold)>log2(1.1)]
length(fold.sig)
significant <- intersect(names(pv.sig), names(fold.sig))
significant.table <- cbind(pv.sig[significant], fold.sig[significant])
significant.table[order(significant.table[,2]),]

# transform pvs
p.trans <- -1 * log10(pv)

# volcano plot
x.line <- -log10(.0001)	#p-value=0.0001
y.line <- log2(1.1)	#fold change=1.1
plot(range(p.trans),range(fold),type='n',xlab=expression(paste("-",log[10]," (p-value)")),ylab=expression(paste(log[2]," fold change")),main='Volcano Plot\nHodgkins Lymphoma')
points(p.trans,fold,col='black',pch=16)
points(p.trans[(p.trans>x.line&fold>y.line)],fold[(p.trans>x.line&fold>y.line)],col=1,pch=21,bg='red')
points(p.trans[(p.trans>x.line&fold<(-1*y.line))],fold[(p.trans>x.line&fold< (-1*y.line))],col=1,pch=21,bg='green')
abline(v=x.line)
abline(h=y.line)
abline(h=(-1* y.line))

#subset data by significant genes
data.sig <- data[significant,]
data.sig.linear <- data.linear[significant,]
dim(data.sig)
dim(data.sig.linear)

#hierarchical clustering, f = failure, s = success of treatment
#TODO: set rownames to f and s
data.transpose <- t(data.sig.linear)
data.transpose.withclass <- data.transpose
rownames(data.transpose.withclass)[gc] <- "s"
rownames(data.transpose.withclass)[act] <- "f"
data.hca <- hclust(dist(data.transpose.withclass,"man"),method="median")
par(mfrow=c(1,1))
plot(data.hca,main='HCA of Hodgkins Lympoma Data \n36 significant genes')

#plot heatmap
heatmap3(as.matrix(t(data.transpose.withclass)),main='2D HCA Hodgkins Lympoma; 36 significant genes',cex.main=0.5)


#PCA for Success vs. Failure
rownames(data.transpose.withclass)[gc] <- "s"
rownames(data.transpose.withclass)[act] <- "f"
data.pca <- prcomp(data.transpose)
data.loads <- data.pca$x[,1:2]
plot(range(data.loads[,1]),range(data.loads[,2]),type="n",xlab='p 1',ylab='p2',main='PCA plot of Lymphoma Data\nTreatment Success vs. Failure') 
points(data.loads[,1][ann=="s"], data.loads[,2][ann=="s"],col=1,bg='red',pch=21,cex=1.5) 
points(data.loads[,1][ann=="f"], data.loads[,2][ann=="f"],col=1,bg='blue',pch=21,cex=1.5) 
legend(0,-2, c("success", "failure"), col=c("red", "blue"), pch=20) 

#PCA for Disease stage: early (1 or 2), late (3 or 4)
metadata@dataTable
class.stage <- data.frame(metadata@dataTable@columns[4]) 
class.stage <- class.stage[-cols.removed,,drop=FALSE]
rownames(class.stage) <- NULL
stage1 <- grep("1", class.stage[,1])
stage2 <- grep("2", class.stage[,1])
stage3 <- grep("3", class.stage[,1])
stage4 <- grep("r4", class.stage[,1])
early <- c(stage1, stage2)
late <- c(stage3, stage4)
length(early)
length(late)
rownames(data.transpose.withclass)[early] <- "e"
rownames(data.transpose.withclass)[late] <- "l"
ann.stage <- rownames(data.transpose.withclass)
plot(range(data.loads[,1]),range(data.loads[,2]),type="n",xlab='p 1',ylab='p2',main='PCA plot of Lymphoma Data\nDisease Stage Early vs. Late') 
points(data.loads[,1][ann.stage=="e"], data.loads[,2][ann.stage=="e"],col=1,bg='red',pch=21,cex=1.5) 
points(data.loads[,1][ann.stage=="l"], data.loads[,2][ann.stage=="l"],col=1,bg='blue',pch=21,cex=1.5) 
legend(0,-2, c("early", "late"), col=c("red", "blue"), pch=20) 

#PCA for gender
metadata@dataTable
class.gender <- data.frame(metadata@dataTable@columns[2]) 
class.gender <- class.gender[-cols.removed,,drop=FALSE]
rownames(class.gender) <- NULL
male <- grep("^male", class.gender[,1])
female <- grep("female", class.gender[,1])
length(male)
length(female)
rownames(data.transpose.withclass)[male] <- "m"
rownames(data.transpose.withclass)[female] <- "f"
data.transpose.withclass[,1:4]
ann.gender <- rownames(data.transpose.withclass)
plot(range(data.loads[,1]),range(data.loads[,2]),type="n",xlab='p 1',ylab='p2',main='PCA plot of Lymphoma Data\nGender: Male vs. Female') 
points(data.loads[,1][ann.gender=="m"], data.loads[,2][ann.gender=="m"],col=1,bg='red',pch=21,cex=1.5) 
points(data.loads[,1][ann.gender=="f"], data.loads[,2][ann.gender=="f"],col=1,bg='blue',pch=21,cex=1.5) 
legend(0,-2, c("male", "female"), col=c("red", "blue"), pch=20) 


# Classification
library(MASS)
data.transpose <- data.frame(ann, data.transpose)
data.transpose
data.sorted <- data.transpose[order(data.transpose[,1]),] #order by class
data.sorted[,1:4]
length(which(data.sorted$ann == "s"))
length(which(data.sorted$ann == "f"))
data.training <- rbind(data.sorted[1:21,], data.sorted[37:92,]) #60%
data.test <- rbind(data.sorted[22:36,], data.sorted[93:127,])
data.training.cl <- as.vector(data.training[,1]) 
data.test.cl <- as.vector(data.test[,1])

# Train on all genes and test
data.lda.all <- lda(data.training.cl~.,data.training[,2:37]) 
data.lda.all 
data.predict.all <- predict(data.lda.all,data.test[,2:37]) 
data.predict.all                          
table(data.predict.all$class,data.test.cl)
 
# plot this
colors = as.numeric(factor(data.predict.all$class)) 
plot(data.predict.all$x,col=colors, pch=21,
     ylab="Discriminant function",axes=F,xlab="Score",      
     main="Discriminant function for Hodgkins Lymphoma Dataset - All Genes") 
axis(1) 
axis(2) 
legend(-0.5,3.0, c("failure", "success"), col=c(1,2), pch=21)
