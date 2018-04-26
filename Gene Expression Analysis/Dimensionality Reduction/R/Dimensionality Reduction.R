data = read.table("c:\\temp\\datasets\\Sotiriou.txt", header=T, row.names=1)
data[1:4,]

data.ann = read.table("C:\\TEMP\\datasets\\Sotiriou_annotations.txt", header=T)
data.ann[1:4,]

dim(data)
colnames(data)

# Cacluate PCA
help(prcomp)
data.pca <- prcomp(t(data))

# plot PCA plot
data.loadings <-data.pca$x[,1:2]
unique(data.ann$site)
plot(range(data.loadings[,1]),range(data.loadings[,2]),type="n",xlab='p1',ylab='p2',main='PCA plot of Sotiriou Data\nKIU vs. OXF sites')
points(data.loadings[,1][data.ann$site=='KIU'], data.loadings[,2][data.ann$site=='KIU'],col=1,bg='red',pch=21,cex=1.5)
points(data.loadings[,1][data.ann$site=='OXF'], data.loadings[,2][data.ann$site=='OXF'],col=1,bg='blue',pch=21,cex=1.5)
legend(0,-10, c("KIU", "OXF"), col=c("red", "blue"), pch=21)

#scree plot
data.pca.var <- round(data.pca$sdev^2 / sum(data.pca$sdev^2)*100,2)
plot(c(1:length(data.pca.var)),data.pca.var,type="b",xlab="# components",ylab="% variance",pch=21,col=1,bg=3,cex=1.5)
title("Scree plot showing % variability explained by each eigenvalue\nKIU/OXF dataset")


#MDS
library(MASS);	
library(multtest);

data.dist <- dist(t(data))

# classical metric MDS on samples (no stress value provided)
data.loc <- cmdscale(data.dist)
plot(data.loc, type="n")
points(data.loc[,1][data.ann$site=='KIU'], data.loc[,2][data.ann$site=='KIU'],bg="red",pch=21,cex=1.5)
points(data.loc[,1][data.ann$site=='OXF'], data.loc[,2][data.ann$site=='OXF'],bg="blue",pch=21,cex=1.5)
title(main='Classic MDS plot of Sotiriou Data\nBy Site')
legend(5,20,c('KIU','OXF'),col=c("red","blue"),pch=19)

# Kruskal’s non-metric MDS on samples
data.mds <- isoMDS(data.dist)
plot(data.mds$points, type = "n")
points(data.mds$points[,1][data.ann$site=='KIU'], data.mds$points[,2][data.ann$site=='KIU'],bg="red",pch=21,cex=1.5)
points(data.mds$points[,1][data.ann$site=='OXF'], data.mds$points[,2][data.ann$site=='OXF'],bg="blue",pch=21,cex=1.5)
title(main='Nonmetric MDS plot of Sotiriou Data\n By Site')
legend(0,20,c('KIU','OXF'),col=c("red","blue"),pch=19)

#center and scale matrix
temp <- t(data)
temp <- scale(temp,center=T,scale=T) 
dim(temp)

# The weighted graph Laplacian
k.speClust2 <- function (X, qnt=NULL) {
  dist2full <- function(dis) {
    n <- attr(dis, "Size")
    full <- matrix(0, n, n)
    full[lower.tri(full)] <- dis
    full + t(full)
  }
  dat.dis <- dist(t(X),"euc")^2
  if(!is.null(qnt)) {eps <- as.numeric(quantile(dat.dis,qnt))}
  if(is.null(qnt)) {eps <- min(dat.dis[dat.dis!=0])}
  kernal <- exp(-1 * dat.dis/(eps))
  K1 <- dist2full(kernal)
  diag(K1) <- 0
  D = matrix(0,ncol=ncol(K1),nrow=ncol(K1))
  tmpe <- apply(K1,1,sum)
  tmpe[tmpe>0] <- 1/sqrt(tmpe[tmpe>0])
  tmpe[tmpe<0] <- 0
  diag(D) <- tmpe
  L <- D%*% K1 %*% D
  X <- svd(L)$u
  Y <- X / sqrt(apply(X^2,1,sum))
}

phi <- k.speClust2(t(temp),qnt=NULL)
plot(range(phi[,1]),range(phi[,2]),xlab="phi1",ylab="phi2",main="Weighted Graph Laplacian plot of Sotiriou Data\nepsilon=0.005")
points(phi[,1][data.ann$site=='KIU'],phi[,2][data.ann$site=='KIU'],col='red',pch=16,cex=1.5)
points(phi[,1][data.ann$site=='OXF'],phi[,2][data.ann$site=='OXF'],col='blue',pch=16,cex=1.5)
legend(-.13,0,c("KIU", "OXF"),col=c('red', 'blue'),pch=15,cex=.7,horiz=F)

