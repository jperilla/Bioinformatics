
#read in classes
classes = read.table("c:\\temp\\eisenClasses.txt", header=T)

#subset the data
dimnames(data)[[2]] #show column names
class1 <- as.character(classes[1:19,2])
data.class1 <- data[,class1]
dimnames(data.class1)[[2]] #show re-ordered column names
class2 <- as.character(classes[20:39,2])
data.class2 <- data[,class2]
