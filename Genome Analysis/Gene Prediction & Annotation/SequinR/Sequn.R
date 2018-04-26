install.packages("seqinr", dependencies = TRUE)
library(seqinr)

# list files in working dir
dir()

# read in fasta doc
palky <- read.fasta(file = "palkyprt.fasta")
palk <- palky[[1]]

# list length, summary table and gc content
length(palk)
table(palk)
GC(palk)

# count dinucleotides
count(palk, 2)

#try with another file
smalty <- read.fasta(file = "smalt.fasta")
smalt <- smalty[[1]]
length(smalt)
table(smalt)
GC(smalt)
count(smalt, 2)
