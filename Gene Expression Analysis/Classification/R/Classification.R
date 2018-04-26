#1
data = read.table("c:\\temp\\datasets\\lung_cancer_data.txt", header=T, row.names=1)
data[,1:4]
colnames(data)

#2
# discriminant analysis
library(MASS)
class <- names(data)
help(rep)
class[grep("Adeno",class)] <- rep("Adeno",length(class[grep("Adeno",class)]))
class[grep("SCLC",class)] <- rep("SCLC",length(class[grep("SCLC",class)]))
class[grep("Normal",class)] <- rep("Normal",length(class[grep("Normal",class)]))
data.transpose <- t(data)
help(data.frame)
data <- data.frame(class, data.transpose)
dim(data)

#3
data.training <- rbind(data[1:6,], data[11:16, ], data[20:22,])
rownames(data.training)
data.test <- rbind(data[7:10,], data[17:19,], data[23:24, ])
rownames(data.test)
data.training.cl <- as.vector(data.training[,1])
data.test.cl <- as.vector(data.test[,1])
data.test <- data.test[,2:3014]

#4
data.lda <- lda(data.training.cl~.,data.training[,2:3])
data.lda
data.predict <- predict(data.lda,data.test[,1:2])
data.predict                         
table(data.predict$class,data.test.cl)

#5
colors = as.numeric(factor(data.predict$class))
plot(data.predict$x,col=colors,
     pch=21,ylab="Discriminant function",axes=F,xlab="Score",
     main="Discriminant function for Lung Cancer Dataset - 2 Genes")
axis(1)
axis(2)
legend(0.5,2.0, c("Adeno", "SCLC", "Normal"), col=c(1,2,3), pch=21)

#6
dim(data.training)
data.lda.ALL <- lda(data.training.cl~.,data.training[,2:3014])
data.predict.ALL <- predict(data.lda.ALL,data.test)
table(data.predict.ALL$class,data.test.cl)

#7

colors = as.numeric(factor(data.predict.ALL$class))
plot(data.predict.ALL$x,col=colors,
     pch=21,ylab="Discriminant function",axes=F,xlab="Score",
     main="Discriminant function for Lung Cancer Dataset - ALL GENES")
axis(1)
axis(2)
legend(0.5,1.0, c("Adeno", "SCLC", "Normal"), col=c(1,2,3), pch=21)
