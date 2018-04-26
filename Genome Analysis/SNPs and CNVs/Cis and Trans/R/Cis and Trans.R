help("read.table")
map = read.table("C:\\TEMP\\datasets\\bp_ct_83_subjects.map")
map
ped = read.table("C:\\TEMP\\datasets\\bp_ct_83_subjects.ped")
ped[1:4,1:10]
data = read.table("C:\\TEMP\\datasets\\bp_ct_u133A.txt")
data

#2
probe <- "216280_s_at"	
snps <- c("SNP_A-2240938","SNP_A-2249234","SNP_A-2120212","SNP_A-1864785","SNP_A-1944204",
          "SNP_A-2101022","SNP_A-2192181","SNP_A-4296300","SNP_A-1961028","SNP_A-2215249",
          "SNP_A-2172180","SNP_A-1945117","SNP_A-1941491")

dicer1 <- data[probe,]
dicer1

dicer.snps <- map[map$V2 %in% snps,]
dicer.snps

dicer.snps.pos <- rownames(dicer.snps)
dicer.snps.pos

#3
m <- matrix(0, ncol = 0, nrow = dim(ped)[[1]])
dicer.snps.ped <- data.frame(m)
for (map.pos in dicer.snps.pos) {
  ped.pos2 <- as.numeric(map.pos)*2+6
  ped.pos1 <- ped.pos2 -1
  colname <- paste("V",as.character(ped.pos1),sep="")
  dicer.snps.ped <- cbind(dicer.snps.ped, ped[,colname])
}

dicer.snps.ped

colnames(dicer.snps.ped) <- dicer.snps$V2
dim(dicer.snps.ped)
dicer.snps.ped

recode <- function(x) {
  x <- as.numeric(x)
  x1 <- seq(1,length(x),by=2)
  x2 <- seq(2,length(x),by=2)
  geno <- paste(x[x1],x[x2],sep="")
  geno[geno=="00"] <- NA
  geno[geno=="11"] <- 0
  geno[geno=="12" | geno=="21"] <- 1
  geno[geno=="22"] <- 2
  geno
}

dim(dicer.snps.ped)
dicer.snps.ped.recoded <- apply(dicer.snps.ped,2,recode)
dicer.snps.ped.recoded

#4
dicer1
dim(dicer1)
dim(dicer.snps.ped.recoded)

lin.mod <- function(x,y) {
  x <- as.numeric(x)
  dat <- data.frame(y,x)
  out <- lm(y~x,data=dat)
  outx <- summary(out)
  return(outx$coefficients["x",4])
}

pvs <- apply(dicer.snps.ped,2,lin.mod,y=as.numeric(dicer1))
pvs

#5

dicer.snps.trans <- map[-(map$V2 %in% snps),]
dicer.snps.trans.14 <- dicer.snps.trans[dicer.snps.trans$V1 == 14,]
dicer.snps.trans.14

dicer.snps.pos.trans <- rownames(dicer.snps.trans.14)
dicer.snps.pos.trans

m <- matrix(0, ncol = 0, nrow = dim(ped)[[1]])
dicer.snps.trans.ped <- data.frame(m)
for (map.pos in dicer.snps.pos.trans) {
  ped.pos2 <- as.numeric(map.pos)*2+6
  ped.pos1 <- ped.pos2 -1
  colname <- paste("V",as.character(ped.pos1),sep="")
  dicer.snps.trans.ped <- cbind(dicer.snps.trans.ped, ped[,colname])
}

dim(dicer.snps.trans.ped)
colnames(dicer.snps.trans.ped) <- dicer.snps.trans.14$V2
dicer.snps.trans.ped

dicer.snps.trans.ped.recoded <- apply(dicer.snps.trans.ped,2,recode)
dicer.snps.trans.ped.recoded

pvs.trans <- apply(dicer.snps.trans.ped,2,lin.mod,y=as.numeric(dicer1))
pvs.trans

#7
help("match")
dicer.snps.ped
cis.dist <- map[map$V2 %in% colnames(dicer.snps.ped),]$V4
cis.dist
cis.pvalue <- -log10(pvs)

trans.dist <- map[map$V2 %in% colnames(dicer.snps.trans.ped),]$V4
trans.dist
trans.pvalue <- -log10(pvs.trans)

par(mfrow=c(2,1)) 
plot(cis.dist, cis.pvalue, main="Cis-acting SNPs log10 p-value vs. distance",
     xlab="Distance", ylab="-log10 p-value")
plot(trans.dist, trans.pvalue, main="Trans-acting SNPs log10 p-value vs. distance",
     xlab="Distance", ylab="-log10 p-value")
