#read in mds data from plink
mds = read.table("C:\\TEMP\\datasets\\hw2\\plink.mds", header=T)
colnames(mds)
mds
#read in clinical data
data.clinical <- read.delim("C:\\temp\\datasets\\hw2\\clinical_table.txt", header=T)
diseases <- data.clinical[,"Profile"]

#plot mds
plot.df <- data.frame(pc1=mds$C1, pc2=mds$C2) 
plot(plot.df, col=as.numeric(diseases), xlab="p1",         
     ylab="p2", main="MDS of disease profiles")
legend(0.03, -0.01, unique(diseases),pch=15,col=unique(as.numeric(diseases)))

#replot MDS by lifetime drug 
drug_use <- data.clinical[,"Lifetime_Drug_Use"]
drug_use

plot(plot.df, col=as.numeric(drug_use), xlab="p1",         
     ylab="p2", main="MDS grouping by Lifetime Drug Use")
legend(0.02, 0.01, unique(drug_use),pch=15,col=unique(as.numeric(drug_use)))

#replot by over or under 50
age <- data.clinical[,"Round_Age"]
age[age < 50] <- 3
age[age >= 50] <- 4
plot(plot.df, col=age, xlab="p1",         
     ylab="p2", main="MDS grouping by Age (over or under 50)")
legend(0.035, -0.02, c(" >= 50", "< 50"),pch=15,col=unique(as.numeric(age)))

left <- data.clinical[,"Left_Brain"]
left

plot(plot.df, col=as.numeric(left), xlab="p1",         
     ylab="p2", main="MDS grouping by Left Brain (Fixed or Frozen)")
legend(0.035, -0.02, unique(left),pch=15,col=unique(as.numeric(left)))

#Lifetime Alcohol 
alcohol_use <- data.clinical[,"Lifetime_Alcohol_Use"]
alcohol_use

plot(plot.df, col=as.numeric(alcohol_use), xlab="p1",         
     ylab="p2", main="MDS grouping by Lifetime Alcohol Use")
legend(0.015, 0.05, unique(alcohol_use),pch=15,col=unique(as.numeric(alcohol_use)))

colnames(data.clinical)

# Loop to run MDS on each parameter
pdf("MDS_Plots.pdf")
for(factor in colnames(data.clinical)) {
  groups <- data.clinical[,factor]
  
  plot(plot.df, col=as.numeric(groups), xlab="p1",         
       ylab="p2", main=paste("MDS grouping by", factor))
  legend(0.01, 0.04, unique(groups),pch=15,col=unique(as.numeric(groups)))
}
dev.off()

plot(plot.df, col=as.numeric(groups), xlab="p1",         
     ylab="p2", main=paste("MDS grouping by", factor))
legend(0.01, -0.015, unique(groups),pch=15,col=unique(as.numeric(groups)))

#4 linear regression analysis of bipolar disorder

# --covar file
mycov <- cbind(as.character(mds$FID),as.character(mds$IID),mds$C1)
mycov <- cbind(mycov, 
               as.integer(data.clinical$sex)-1, 
               as.integer(data.clinical$Left_Brain)-1)
mycov
write.table(mycov,"mycov.txt",sep=" ",row.names=F,col.names=F,quote=F)

# --keep file
data.clinical$Database_ID
mds$FID
data.clinical$Profile
bip_control <- data.clinical[data.clinical$Profile=="BP" |
                    data.clinical$Profile=="Unnaffected control",]$Database_ID
bip_control <- cbind(as.character(bip_control),as.character(mds[mds$FID %in% bip_control,]$IID))
bip_control
write.table(bip_control,"bipolar_control.txt",sep=" ",row.names=F,col.names=F,quote=F)

#--pheno file (1=unaffected, 2=affected)
pheno <- cbind(as.character(mds$FID),mds$IID)
temp <- data.clinical
temp$Profile
temp$Profile <- as.integer(temp$Profile)
temp[!temp$Profile==4,] <- 2
temp[temp$Profile==4,] <- 1
pheno <- cbind(pheno, temp$Profile)
pheno
write.table(pheno,"pheno.txt",sep=" ",row.names=F,col.names=F,quote=F)

# --keep file for schiz
data.clinical$Profile
schz_control <- data.clinical[data.clinical$Profile=="Schiz." |
                               data.clinical$Profile=="Unnaffected control",]$Database_ID
schz_control <- cbind(as.character(schz_control),as.character(mds[mds$FID %in% schz_control,]$IID))
schz_control
write.table(schz_control,"schz_control.txt",sep=" ",row.names=F,col.names=F,quote=F)

#read in assoc
bip_summary <- read.table("C:\\TEMP\\datasets\\hw2\\plink.assoc.logistic", header=T)
plessthan01.bip <- bip_summary[which(bip_summary$TEST=="ADD" & bip_summary$P<.01),]
plessthan01.bip

schiz_summary <- read.table("C:\\TEMP\\datasets\\hw2\\plink.assoc.logistic.schiz", header=T)
plessthan01.schiz <- schiz_summary[which(schiz_summary$TEST=="ADD" & schiz_summary$P<.01),]
plessthan01.schiz
