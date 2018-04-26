# read in data
data = read.table("c:\\temp\\spellman.txt", header=T, row.names=1)

#show subset of data
data[1:4, 1:4]

# show dimensions of data
dim(data)

# show names of samples
dimnames(data)[[2]]

# subset columns
cdc15Data = data[,23:46]

# subset rows
cdc15Data = data[1:6,]
