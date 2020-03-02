##
##  This code computes for differentially expressed genes using DESeq2
##

#Import package DESeq2
library("DESeq2")
# Create the DESeqDataSet object from Matrix of counts and metadata
dds <- DESeqDataSetFromMatrix(countData = round(df),colData = metadatah,design = ~Groups)
nrow(dds) 

#Run DESeq function on the data to perform 
#differential gene expression analysis
dds4 <- DESeq(dds)
colData(dds4)

# Building out results table
res <- results(dds4)
mcols(res , use.names=TRUE)
summary(res )

#Orderin te p values accordin to sinificance
res11  <- res [order(res$pvalue),]
sum(res11$padj < 0.1, na.rm=TRUE)

#We order our results table by the smallest p value:
resSig <- subset(res, padj < 0.1)
nrow(resSig)

write.table(as.data.frame(rownames(resSig)),paste0("deseq2",j,".csv"),row.names = F,sep=",",col.names=FALSE)
df1=read.csv(paste0("deseq2",j,".csv"),header = F)
source('putinalist.R')
putinalist(df1)
