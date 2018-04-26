#3 - Read in filtered data
data = read.table("~//Datasets//FC305JN_cols.txt")
dim(data)

#filter data for MACs
colnames(data) <- c("Name", "Sequence", "chr6.fa", "Locus", "Strand")
colnames(data)
data.F <- data[data$Strand=="F",]
numUnique.F <- length(unique(data.F$Locus))
data.R <- data[data$Strand=="R",]
numUnique.R <- length(unique(data.R$Locus))
numDuplicates <- dim(data)[[1]] - (numUnique.F + numUnique.R)

#output MACS file
macs <- data.frame("Chr"=as.character(), "start"=as.numeric(), "end"=as.numeric(), 
                   "??"=as.numeric(), "tags"=as.numeric(), "sense"=as.character())

for (row in 1:nrow(data)) {
  print(row)
  start <- data[row,4]
  if(nrow(macs[macs$start==start,]) > 0) {
    macs[macs$start==start,5] <- macs[macs$start==start,5] + 1
  } else {
    end <- as.numeric(data[row, 4]) + nchar(as.character(data[row, 2]))
    strand <- if(data[row, 5]=="F") "+" else "-"
    new <- data.frame("chr6", start, end, 0, 1, strand)
    names(new) <- c("Chr", "start", "end", "??", "tags", "sense")
    macs <- rbind(macs, new)
  }
}

macs.sorted <- macs[order(macs$start),]
write.table(macs.sorted, file="~//Datasets//FC305JN_macs.txt",
            row.names=FALSE, quote=FALSE, sep='\t', col.names=FALSE)

#4
macs.data <- read.table("~//Datasets//ame_peaks.xls", header=T)
summary(macs.data)

plot(macs.data$tags, macs.data$start, type="n",
     xlab="", ylab="No. of Tags", xaxt="n",
     ylim=c(0,100),
     xlim=c(0,max(macs.data$start)),
     main="Chip-seq peaks found using MACS\nChromosome 6", 
     col = "black")

lines(macs.data$start, macs.data$tags, type="h", col="red") 
lines(macs.sorted$start, macs.sorted$tags, type="h", col="blue")

minStart = min(macs.data$start)
maxStart = max(macs.data$start)
range = maxStart - minStart
xticks<-seq(minStart, maxStart, by=round(range/20)) 
axis(side=1, at=xticks, labels=xticks,las=2, cex.axis=0.5)
mtext("Position", side=1, line=4)
