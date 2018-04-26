
data = read.table("c:\\temp\\datasets\\agingStudy11FCortexAffy.txt", header=T, row.names=1)
data[1:4,]

data.ann = read.table("C:\\TEMP\\datasets\\agingStudy1FCortexAffyAnn.txt", header=T)
data.ann[]

#2
g.g <- c(1394,  1474,  1917,  2099,  2367,  2428, 2625,  3168,  3181,  3641,  3832,  4526,
         4731,  4863,  6062,  6356,  6684,  6787,  6900,  7223,  7244,  7299,  8086,  8652,
         8959,  9073,  9145,  9389, 10219, 11238, 11669, 11674, 11793)
data.gender <- data[g.g,]
dim(data.gender)
data.gender.male <- data.gender[,1:18]
data.gender.female <- data.gender[,19:30]

g.a <- c(25, 302,  1847,  2324,  246,  2757, 3222, 3675,  4429,  4430,  4912,  5640, 5835, 5856,  6803,  7229,  7833,  8133, 8579,  8822,  8994, 10101, 11433, 12039, 12353,
         12404, 12442, 67, 88, 100)
data.age <- data[g.a,]

data.ann <- data.ann[order(data.ann$Age),]
data.ann.50andup <- as.vector(data.ann[data.ann$Age >= 50, 1])
data.ann.under50 <- as.vector(data.ann[data.ann$Age < 50, 1])

colnames(data.age) <- lapply(colnames(data.age), 
                              function(x) { strsplit(x, ".", fixed=TRUE)[[1]][1]})
data.age.50andup <- data.age[,data.ann.50andup]
data.age.under50 <- data.age[,data.ann.under50]
data.age.sorted <- cbind(data.age.under50, data.age.50andup)


#3
t.test.all.genes <- function(x,s1,s2) {
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- t.test(x1,x2, alternative='two.sided',var.equal=T)
  out <- as.numeric(t.out$p.value)
  return(out)
}

rawp.gender <- apply(data.gender,1, t.test.all.genes,s1=1:18,s2=19:30)
rawp.age <- apply(data.age.sorted, 1, t.test.all.genes, s1=1:12, s2=13:30)

source("https://bioconductor.org/biocLite.R")
biocLite("multtest")
library(multtest) 

help(mt.rawp2adjp)
help(mt.plot)
adjustedp.gender <- mt.rawp2adjp(rawp.gender, proc="Holm")
adjustedp.age <- mt.rawp2adjp(rawp.age, proc="Holm")

#4
mt.plot(adjustedp.gender$adjp,plottype="pvsr",proc=c("Holm"),
        lwd=2, leg=c(1,0.08))
title("Male and Female\n Holm's adjusted p-value \nvs. number of rejected hypotheses") 
adjustedp.age$adjp
mt.plot(adjustedp.age$adjp,plottype="pvsr",proc=c("Holm"),
        lwd=2, leg=c(1,0.8))
title("Under 50 and 50 and up \nHolm adjusted p-value \nvs. number of rejected hypotheses")

#6
bonf.gender <- mt.rawp2adjp(rawp.gender, proc="Bonferroni")
bonf.age <- mt.rawp2adjp(rawp.age, proc="Bonferroni")
bonf.gender$adjp
mt.plot(bonf.gender$adjp,plottype="pvsr",proc=c("Bonferroni"),
        lwd=2, leg=c(1,0.08))
title("Male and Female\n Bonferroni adjusted p-value \nvs. number of rejected hypotheses") 
mt.plot(bonf.age$adjp,plottype="pvsr",proc=c("Bonferroni"),
        lwd=2, leg=c(1,0.8))
title("Under 50 and 50 and up \nBonferroni adjusted p-value \nvs. number of rejected hypotheses")


#Part 2
data = read.table("c:\\temp\\datasets\\tcga_brca_fpkm.txt", header=T, row.names=1)
data[1:4,1:4]
dim(data)

help(read.table)
data.ann = read.table("C:\\TEMP\\datasets\\tcga_brca_fpkm_sam.txt", header=T, fill=T)
data.ann

help(grep)
grep("^*GATA3&?", rownames(data), ignore.case=FALSE)
data.gata3 <- as.numeric(data[6362,])
data.gata3
length(data.gata3)
#get the number that each need to be, to be in the top 25% 
twentyfive <- round(119/4)
topthreshold <- sort(data.gata3, decreasing=T)[twentyfive]
group <- lapply(data.gata3, function(x) { if(x > topthreshold) { return (1);} else {return (0);} })
group <- as.numeric(group)


data.ann.withgroup <- cbind(data.ann, group)
data.ann.withgroup

#11
set_status <- function(x) {
  if (is.na(x)) {return (NA)}
  
  if (x == "DECEASED") {return (1)}
  else {return (0)}
}

status <- unlist(lapply(data.ann.withgroup[,"vital_status"],set_status))
time <- as.vector(data.ann.withgroup[,"months_to_event"])
x <- data.frame(status,time)
names(x) <- c("time","status")
library(survival)
survdiff(Surv(time, status) ~ group)

#12
help(coxph.object)
fit <- coxph(Surv(time, status) ~ group)
summary(fit)

#13
f<-survfit(Surv(time, status) ~ group,type="kaplan-meier")
colors = c("blue", "green")
plot(f, lty = 2:3, xlab="Months", ylab="S(t)", col=colors)
legend(90, .9, c("Lower 75 percentile", "Top 25 percentile"), lty = 2:3, col=colors)
legend(5,0.3, c("KM p = 0.0788", "Cox p=0.06188", "Hazard Ratio = 2.309"), )
title("Kaplan-Meier Curves\nfor Top 25 Percentile Expression")


