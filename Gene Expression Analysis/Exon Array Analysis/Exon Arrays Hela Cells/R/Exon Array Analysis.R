# 2 - read in data
e = read.table("c:\\temp\\datasets\\homework\\exon-rma-sketch.summary.txt", header=T, row.names=1)
dim(e)

g = read.table("c:\\temp\\datasets\\homework\\gene-rma-sketch.summary.txt", header=T, row.names=1)
dim(g)

p = read.table("c:\\temp\\datasets\\homework\\dabg.summary.txt", header=T, row.names=1)
dim(p)

map = read.csv("c:\\temp\\datasets\\homework\\HuEx-1_0-st-v2.na24.hg18.probeset_abbr.csv", header=T, row.names=1)
dim(map)
map

#3 - PCA to find outliers
e.pca <- prcomp(t(e))
e.loads <- e.pca$x[,1:2]
plot(e.loads[,1],e.loads[,2],main="2D PCA plot - Exons",
     xlab="p1",ylab="p2",col='red',cex=1,pch=15)
text(e.loads,label=dimnames(e)[[2]],pos=3,cex=0.5)


g.pca <- prcomp(t(g))
g.loads <- g.pca$x[,1:2]
plot(g.loads[,1],g.loads[,2],main="2D PCA plot - Genes",
     xlab="p1",ylab="p2",col='red',cex=1,pch=15)
text(g.loads,label=dimnames(g)[[2]],pos=3,cex=0.5)

#4
#intersect matching probes between annotation file and data matrix
x <- intersect(dimnames(e)[[1]],dimnames(map)[[1]])
length(x)

#subset the rows in the exon, p-value matrix, and annotation file to the intersecting probes
e.dat <- e[x,]
dim(e.dat)
p.dat <- p[x,]
dim(p.dat)
map.dat <- map[x,]
dim(map.dat)

e.dat[1:4,]

#get unique transcript cluster IDs (gene IDs) from annotation file
u <- unique(as.character(map.dat$transcript_cluster_id))
u <- intersect(u,dimnames(g)[[1]])
length(u)

# two t-test called by exon ni function below
t.two <- function(x,sam,v=F) {
  x <- as.numeric(x)
  out <- t.test(as.numeric(x[sam]),as.numeric(x[!sam]),alternative="two.sided",var.equal=v)
  control <- mean(log2(x[sam]), na.rm=TRUE)
  test <- mean(log2(x[!sam]), na.rm=TRUE)
  fold <- control - test
  print(fold)
  o <- as.numeric(c(out$statistic,out$p.value,out$conf.int[1],out$conf.int[2], fold))
  names(o) <- c("test_statistic","pv","lower_ci","upper_ci", "fold_change")
  return(o)
}

# exon ni function
exon.ni <- function(genex,exonx,rx) {
  ni <- t(t(exonx)-genex)
  ttest <- t(apply(ni,1,t.two,sam=as.logical(rx),v=F))
  return(ttest)	
}

# treated samples vs controls
pvalues.clotrim <- data.frame(test_statistic=as.numeric(), pv=as.numeric(),lower_ci=as.numeric(), upper_ci=as.numeric(), 
                              fold_change=as.numeric())
r.clotrim = c(0,0,0,1,1,1)

pvalues.chlorex <- data.frame(test_statistic=as.numeric(), pv=as.numeric(),lower_ci=as.numeric(), upper_ci=as.numeric(), 
                              fold_change=as.numeric())
r.chlorex = c(0,0,0,1,1,1)
count <- 1
for (uniqId in u) {
  print(count)
  ex <- dimnames(map.dat[map.dat$transcript_cluster_id %in% uniqId,])[[1]]
  d.exon.clotrim <- e.dat[ex,1:6]
  d.gene.clotrim <- g[uniqId,1:6]	
  if(dim(d.exon.clotrim)[[1]] > 2) {
    ni.out.clotrim <- exon.ni(genex=as.numeric(d.gene.clotrim),exonx=d.exon.clotrim,rx=r.clotrim)
    pvalues.clotrim <- rbind(pvalues.clotrim, ni.out.clotrim)
  }
  
  d.exon.chlorex <- e.dat[ex,10:15]
  d.gene.chlorex <- g[uniqId,10:15]	
  if(dim(d.exon.chlorex)[[1]] > 2) {
    ni.out.chlorex <- exon.ni(genex=as.numeric(d.gene.chlorex),exonx=d.exon.chlorex,rx=r.chlorex)
    pvalues.chlorex <- rbind(pvalues.chlorex, ni.out.chlorex)
  }
  count <- count + 1
}


dim(pvalues.chlorex)
dim(pvalues.clotrim)
length(u)

pvalues.clotrim


#5 Find significant p-values and |fold change| > 1.5, remember
# to retranspose fold change
clotrim.sigp <- na.omit(pvalues.clotrim[pvalues.clotrim$pv < 0.01,])
dim(clotrim.sigp)

clotrim.sigfc <- na.omit(pvalues.clotrim[abs(2^pvalues.clotrim$fold_change) > 1.5,])
dim(clotrim.sigfc)

clotrim.sig <- intersect(dimnames(clotrim.sigp)[[1]],dimnames(clotrim.sigfc)[[1]])
length(clotrim.sig)

######
chlorex.sigp <- na.omit(pvalues.chlorex[pvalues.chlorex$pv < 0.01,])
dim(chlorex.sigp)

chlorex.sigfc <- na.omit(pvalues.chlorex[abs(2^pvalues.chlorex$fold_change) > 1.5,])
dim(chlorex.sigfc)

chlorex.sig <- intersect(dimnames(chlorex.sigp)[[1]],dimnames(chlorex.sigfc)[[1]])
length(chlorex.sig)

