#1
library(multtest)
library(Biobase)
library(annotate)

data(golub);

#2
data <- data.frame(golub)
data[1:4,]


#3
ann.dat2 <- golub.cl

#4
help("wilcox.test")
wilcox.test.all.genes <- function(x,s1,s2) {
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- wilcox.test(x1,x2, alternative='two.sided', exact=F, correct=T)
  out <- as.numeric(t.out$statistic)
  return(out)
}

original.wmw.run <- apply(data,1,wilcox.test.all.genes,s1=ann.dat2==0,s2=ann.dat2==1)
max(original.wmw.run)

#5
length(colnames(data[ann.dat2==0]))
length(colnames(data[ann.dat2==1]))
maxteststats <- list()
iterate <- for (i in 1:500) {
  x1 <- sample(colnames(data), 27)
  x2 <- sample(colnames(data), 11)
  wmw.run <- apply(data,1,wilcox.test.all.genes,s1=x1,s2=x2)
  maxteststats[i] <- max(wmw.run)
}

iterate()

max(as.numeric(maxteststats))
max(original.wmw.run)
original.wmw.run[1:30]

#6
max <- max(as.numeric(maxteststats))
ninetyfive <- max * .95
top.ninetyfive <- original.wmw.run[original.wmw.run > ninetyfive]
top.ninetyfive

#7
library(limma)
help("ebayes")

design <-cbind(Grp1=1,Grp2vs1=c(rep(1,27),rep(0,11)))
fit <- lmFit(data,design)
fit <- eBayes(fit)
pv <- fit$p.value[,2]
pv

#8
sortedps <- sort(pv)
n <- length(top.ninetyfive)
lowest.n <- sortedps[1:n]
lowest.n
common <- intersect(names(lowest.n), names(top.ninetyfive))
length(common)

#9
t.test.all.genes <- function(x,s1,s2) {
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- t.test(x1,x2, alternative='two.sided', exact=F, correct=T)
  out <- as.numeric(t.out$p.value)
  return(out)
}

students <- apply(data,1,t.test.all.genes,s1=ann.dat2==0,s2=ann.dat2==1)
students.lessthan <- students[students < 0.01]
pv[names(students.lessthan)]
plot(pv[names(students.lessthan)],students.lessthan,
     xlab='empirical Bayes',ylab='Student’s t-test',
     main='P-value distribution comparison',
     cex=0.5,col=3,pch=15)

