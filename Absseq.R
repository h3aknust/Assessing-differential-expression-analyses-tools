
##
#This code computes for differentially expressed genes using ABSSeq
##
#Import ABSSeq
library(ABSSeq)

#create
data(simuN5)
obj <- ABSDataSet(simuN5$counts, factor(simuN5$groups))
obj <- ABSSeq(obj)
res <- results(obj,c("Amean","Bmean","foldChange","pvalue","adj.pvalue"))
#head(res)

data(simuN5)
groups<-factor(simuN5$groups)
obj <- ABSDataSet(counts=simuN5$counts)
design <- model.matrix(~0+groups)
res <- ABSSeqlm(obj,design,condA=c("groups0"),condB=c("groups1"))
###############SET_1################################################################################
#set conditions
group <- as.character(metadatah$Groups)

obj <- ABSDataSet(df, factor(group))
obj <- normalFactors(obj)

head(counts(obj,norm=TRUE))
obj=callParameter(obj)
obj <- callDEs(obj)

ABss <- (results(obj,c("pvalue","adj.pvalue")))
sum(ABss[,2] <  0.05)
ABss_DE1 <- subset(ABss, ABss[,2] <  0.05)
ABss_DE1 <- rownames(ABss_DE1)

write.csv(ABss_DE1, file = "ABss_DE1.csv", row.names = FALSE)
write.table(as.data.frame(ABss_DE1),paste0("ABss",j,".csv"),row.names = F,sep=",",col.names=FALSE)
df1=read.csv(paste0("ABss",j,".csv"),header = F)
source('putinalist.R')
putinalist(df1)

#####Afold
obja <- ABSDataSet(counts=df, groups=factor(group))
obja <- ABSSeq(obja, useaFold=TRUE)
resa <- results(obja,c("Amean","Bmean","foldChange","pvalue","adj.pvalue"))
ABssa <- (results(obja,c("pvalue","adj.pvalue")))

ABssa_DE1 <- subset(ABssa, ABssa[,2] <  0.05)
ABssa_DE1 <- rownames(ABssa_DE1)

#write.csv(ABssa_DE1, file = "ABssa_DE1.csv", row.names = FALSE)
write.table(as.data.frame(ABssa_DE1),paste0("ABssa",j,".csv"),row.names = F,sep=",",col.names=FALSE)
df1=read.csv(paste0("ABssa",j,".csv"),header = F)
source('putinalist.R')
putinalist(df1)

