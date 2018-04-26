# 1 - read in data
e = read.table("c:\\temp\\datasets\\exon-rma-sketch.summary.txt", header=T, row.names=1)
e[1:4,1:4]

g = read.table("c:\\temp\\datasets\\gene-rma-sketch.summary.txt", header=T, row.names=1)
g[1:4,1:4]

p = read.table("c:\\temp\\datasets\\dabg.summary.txt", header=T, row.names=1)
p[1:4,1:4]

map = read.csv("c:\\temp\\datasets\\HuEx-1_0-st-v2.na24.hg18.probeset_abbr.csv", header=T, row.names=1)
map[1:4,]

#2 
#class membership
r <- c(1,0,1,1,1,0,1,1,1,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1)

#intersect matching probes between annotation file and data matrix
x <- intersect(dimnames(e)[[1]],dimnames(map)[[1]])

#subset the rows in the exon, p-value matrix, and annotation file to the intersecting probes
e <- e[x,]
p <- p[x,]
map <- map[x,]

#get unique transcript cluster IDs (gene IDs) from annotation file
u <- unique(as.character(map$transcript_cluster_id))
u <- intersect(u,dimnames(g)[[1]])

#3 filter exon probes based on p-values in dabg (p)
count.det <- function(x) { length(which(x>0.05))}
det <- apply(p, 1, count.det)
probes_to_remove <- names(det[det > 13]) # filter these out
p_filtered <- p[-which(rownames(p) %in% probes_to_remove),]
e_filtered <- e[-which(rownames(e) %in% probes_to_remove),]
dim(e_filtered)
dim(p_filtered)

#4 filter map and intersect with gene table
map_filtered <- map[-which(rownames(map) %in% probes_to_remove),]
u_filt <- unique(as.character(map_filtered$transcript_cluster_id))
u_filt <- intersect(u_filt,dimnames(g)[[1]])
length(u_filt)

# 5
# two t-test called by exon ni function below
t.two <- function(x,sam,v=F) {
  x <- as.numeric(x)
  out <- t.test(as.numeric(x[sam]),as.numeric(x[!sam]),alternative="two.sided",var.equal=v)
  o <- as.numeric(c(out$statistic,out$p.value,out$conf.int[1],out$conf.int[2]))
  names(o) <- c("test_statistic","pv","lower_ci","upper_ci")
  return(o)
}

# exon ni function
exon.ni <- function(genex,exonx,rx) {
  ni <- t(t(exonx)-genex)
  ttest <- t(apply(ni,1,t.two,sam=as.logical(rx),v=F))
  return(ttest)	
}


pvalues <- data.frame(test_statistic=as.numeric(), pv=as.numeric(),
                      lower_ci=as.numeric(), upper_ci=as.numeric())
for (uniqId in u_filt) {
  ex <- dimnames(map[map$transcript_cluster_id %in% uniqId,])[[1]]
  d.exon <- e[ex,]
  d.gene <- g[uniqId,]	
  if(dim(d.exon)[[1]] > 2) {
    ni.out <- exon.ni(genex=as.numeric(d.gene),exonx=d.exon,rx=r)
    pvalues <- rbind(pvalues, ni.out)
  }
}

#6 sort by p-value and get transcript cluster with lowest p-value
lowest <- pvalues[order(pvalues$pv),][1,]
lowest$pv
sig_exon <- rownames(lowest)
sig_transcript <- as.character(map[rownames(map)==sig_exon, ]$transcript_cluster_id)

#7 
#boxplot of exons for gene 2426958 (sig_transcript)
plot.exons <- function(exonx,genex,rx,ti) {
  rr <- rx
  rx <- rep(rx,nrow(exonx))
  rx[rx==1] <- "A"
  rx[rx==0] <- "B"
  rx <- as.factor(rx)
  ni <- t(t(exonx)-genex)
  exonx <- as.data.frame(t(ni))
  ex.stack <- stack(exonx)
  d <- data.frame(ex.stack,rx)
  names(d) <- c("exon_values","exon_id","class")
  
  d$exon_id <- as.factor(d$exon_id)
  d$class <- as.factor(d$class)
  genex.title <- as.character(map[match(ti,as.character(map$transcript_cluster_id)),"gene_assignment"])
  plot(c(.5,(ncol(exonx)+.5)),range(d[,1]),type="n",axes=F,xlab="",ylab="")
  boxplot(exon_values~exon_id,add=T,subset=d$class=="A",d,col="salmon",border='red',cex.axis=.75,las=2,ylab='Log2 normalized intensity',main=paste("Gene ID:",ti,"\n",genex.title),boxwex=0.4)
  boxplot(exon_values~exon_id,subset=d$class=="B",d,add=T,col="green",border='darkgreen',axes=F,boxwex=0.4, at=c(1:ncol(exonx))+0.1)
}

d.exon <- e[sig_exon,]
d.gene <- g[sig_transcript,]	
plot.exons(exonx=d.exon,genex=as.numeric(d.gene),rx=r,ti=sig_transcript)


