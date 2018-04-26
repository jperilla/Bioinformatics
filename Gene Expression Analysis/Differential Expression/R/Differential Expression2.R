
data = read.table("c:\\temp\\rat_KD\\rat_KD.txt", header=T, row.names=1)

# log2 the data
data.log2 <- log2(data)
dim(data.log2)[1]

# function to calculate Student’s two-sample t-test on all genes at once
# function returns the p-value for the test
# NAs are removed for each test
# s1 and s2 are dimensions of the two samples
# run function on each gene in the data frame
t.test.all.genes <- function(x,s1,s2) {
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- t.test(x1,x2, alternative="two.sided",var.equal=T)
  out <- as.numeric(t.out$p.value)
  return(out)
}

names(data.log2)
dim(data.log2)
control <- 1:6
keto <- 7:11
pv <- apply(data.log2,1,t.test.all.genes,s1=control,s2=keto)

#plot a histogram of the p-values
# look at distribution of p-values
par(mfrow=c(1,2))
hist(pv,col="lightblue",xlab="p-values",
     main="P-value dist’n between\ncontrol diet and \nketogenic diet groups",
     cex.main=0.9)
abline(v=.05,col=2,lwd=2)
hist(-log10(pv),col="lightblue",xlab="log10(p-values)", 
     main="-log10(pv) dist’n between\ncontrol diet and \nketogenic diet groups",cex.main=0.9)
abline(v=-log10(.05),col=2,lwd=2)

length(pv)
pv.lessthan05 <- pv[pv < .05]
length(pv.lessthan05)

pv.lessthan01 <- pv[pv < .01]
length(pv.lessthan01)

alpha <- 0.5/15923
alpha
pv.lessthanalpha <- pv[pv < 0.5/15923]
length(pv.lessthanalpha)
sum(pv<(.05/15923))

help(sum)

# calculate means of each class for each gene
control.m <- apply(data.log2[,control],1,mean,na.rm=T)
keto.m <- apply(data.log2[,keto],1,mean,na.rm=T)

# get fold changes by subtracting one from the other
fold <- control.m-keto.m

# max and min fold change, convert back to linear
help(max)
help(exp)
fold.max <- 2^max(fold)
fold.max
fold.min <- 2^min(fold)
fold.min

log2(fold.min)

# find p-values less than threshold, and |fold change| > 2
pv.sig <- pv.lessthanalpha
fold <- 2^fold
fold.sig <- fold[abs(fold)>2]
names(fold.sig)
significant <- intersect(names(fold.sig), names(pv.sig))
significant

# transpose p-values
pv.trans <- -1 * log10(pv)
pv.trans

# transpose fold back to log2
fold <- log2(fold)
max(fold)

# volcano plot pv.trans vs. log2(fold)
plot(range(pv.trans),range(fold),type='n',
     xlab='-1*log10(p-value)',ylab='fold change',
     main='Volcano Plot\nControl and Ketogenic Diet\ngroup differences')
points(pv.trans,fold,col='black',pch=20,bg=1)
points(pv.trans[(pv.trans> -log10(.05)&fold>log2(2))],fold[(pv.trans> -log10(.05)&fold>log2(2))],col=1,bg=2,pch=21)
points(pv.trans[(pv.trans> -log10(.05)&fold< -log2(2))],fold[(pv.trans> -log10(.05)&fold< -log2(2))],col=1,bg=3,pch=21)
abline(v= -log10(.05))
abline(h= -log2(2))
abline(h=log2(2))

