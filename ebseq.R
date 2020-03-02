##
#This code computes for differentially expressed genes using EBSeq
##
#Loading the EBSeq package
library(EBSeq)

#Estimating the library sizes with converted factor to Numeric
librarySizes = MedianNorm(df)


####EBSeq take matrix entry not dataframe... so convert dataframe into matrix
df1=as.matrix(df)

##Defining the conditions
Conditions1=as.factor(rep(c("sample1","sample2"),each=nrow(metadatah)/2))

##computing The DGE with EBTest
EBOut=EBTest(Data=df1,Conditions=Conditions1,sizeFactors=librarySizes, maxround=5)

viewResults = GetPPMat(EBOut)

##Extracting only EBSeq Results
EBDERes=GetDEResults(EBOut)


##Extracting the DGE using a false discovery rate of 0.05
EBDERes=GetDEResults(EBOut, FDR=0.05)

write.table(as.data.frame(EBDERes$DEfound),paste0("ebseq",j,".csv"),row.names = F,sep=",",col.names=FALSE)
df1=read.csv(paste0("ebseq",j,".csv"),header = F)
source('putinalist.R')
putinalist(df1)


