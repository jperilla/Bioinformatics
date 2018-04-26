

#globa median normalization
swirl3.norm <- maNorm(swirl3, norm="median")
plot(swirl3.norm, lines.func=NULL, legend=NULL,
     main="Swirl array 3: post--normalization MA--plot")

help(maNorm)
swirl3.norm.gi <- maNorm(swirl3, norm="loess")
plot(swirl3.norm.gi, lines.func=NULL, legend=NULL,
     main="Swirl array 3: post--global intensity\nnormalization MA--plot")

dir.path <- "C:\\temp\\"
a.cdna <- read.GenePix(path=dir.path,name.Gf = "F532 Median",name.Gb ="B532 Median", name.Rf = "F635 Median", name.Rb = "B635 Median",name.W ="Flags")
a.cdna

a.cdna.patient1 <- a.cdna[,1]
a.cdna.patient2 <- a.cdna[,2]
par(mfrow=c(3,1)) 
plot(a.cdna.patient1,lines.func=NULL, legend=NULL,
     main="Patient 1: pre--normalization MA--plot")

a.cdna.patient1.printeTipLoess <- maNorm(a.cdna.patient1, norm="printTipLoess")
plot(a.cdna.patient1.printeTipLoess, lines.func=NULL, legend=NULL,
     main="Patient1: post-normalization MA--plot\n print Tip Loess")

a.cdna.patient1.scalePrintTipMAD <- maNorm(a.cdna.patient1, norm="scalePrintTipMAD")
plot(a.cdna.patient1.scalePrintTipMAD, lines.func=NULL, legend=NULL,
     main="Patient1: post-normalization MA--plot\n scale Print Tip MAD")


par(mfrow=c(3,1)) 
plot(a.cdna.patient2,lines.func=NULL, legend=NULL,
     main="Patient 2: pre--normalization MA--plot")

a.cdna.patient2.printeTipLoess <- maNorm(a.cdna.patient2, norm="printTipLoess")
plot(a.cdna.patient2.printeTipLoess, lines.func=NULL, legend=NULL,
     main="Patient2: post-normalization MA--plot\n print Tip Loess")

a.cdna.patient2.scalePrintTipMAD <- maNorm(a.cdna.patient2, norm="scalePrintTipMAD")
plot(a.cdna.patient2.scalePrintTipMAD, lines.func=NULL, legend=NULL,
     main="Patient2: post-normalization MA--plot\n scale Print Tip MAD")


#make a data frame to export
help(maNorm)
p1.print <- maM(a.cdna.patient1.printeTipLoess)
p2.print <- maM(a.cdna.patient2.printeTipLoess)
probesets <- maLabels(maGnames(a.cdna))
df.print <- cbind(p1.print, p2.print)
dimnames(df.print)[[1]] <- probesets
df.print

p1.scale <- maM(a.cdna.patient1.scalePrintTipMAD)
p2.scale <- maM(a.cdna.patient2.scalePrintTipMAD)
df.scale <- cbind(p1.scale, p2.scale)
dimnames(df.scale)[[1]] <- probesets
df.scale

source("https://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("limma")

source("https://bioconductor.org/biocLite.R")
biocLite("simpleaffy")
biocLite("affyPLM")
biocLite("fpc")

fns <- sort(list.celfiles(path="c:/temp/",full.names=TRUE))
data.affy <- ReadAffy(filenames=fns,phenoData=NULL)
data.affy

library(simpleaffy)
data.affy.mas <- justMAS(data.affy)
data.affy.mas
mas <- exprs(data.affy.mas)
dim(mas)
data.affy.rma <- rma(data.affy)
data.affy.rma
rma <- exprs(data.affy.rma)
dim(rma)
