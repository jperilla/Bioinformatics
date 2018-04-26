#1
source("https://bioconductor.org/biocLite.R")
biocLite("fibroEset")
library(fibroEset)
data(fibroEset)
data.cl <- fibroEset$species
data <- exprs(fibroEset)
data[1:4,1:4]

#2
genesample <- sample(rownames(data), 50)
genesubset <- data[genesample,]

help(dist)
help(hclust)

#3
colnames(genesubset) <- data.cl
genes.transposed <- t(genesubset)
genes.dist.man <- dist(genes.transposed,method="manhattan")
genes.clust.man <- hclust(genes.dist.man, method="median")
genes.clust.man
plot(genes.clust.man, 
     main="Hierarchical Manhattan clustering of 50 genes - fibroEset")

#4
hm.rg <- c("#FF0000","#CC0000","#990000","#660000","#330000","#000000","#000000","#0A3300","#146600","#1F9900","#29CC00","#33FF00")
heatmap(as.matrix(genesubset),col=hm.rg, 
        main="HCA of 50 genes - fibroEset")

#5
genes.pca <- prcomp(genesubset)
data.loads <- genes.pca$x[,1:2]
cl <- kmeans(data.loads,centers=3)	
cl$centers
length(cl$cluster)
length(genes.pca$rotation[,1])
dim(genes.pca$x)

#6
genes.pca$x[,1:2]

plot.df <- data.frame(pc1=data.loads[,1], 
                      pc2=data.loads[,2])

plot(plot.df, col=cl$cluster, xlab="Eigenvector 1",
        ylab="Eigenvector 2", main="k-means classification using \nPCA eigenvector 1 and eigenvector 2")
