##
#This code computes for differentially expressed genes using bayseq
##
#Import package bayseq
library(baySeq)

# Specify the number of clusters to run
if(require("parallel")) cl <- makeCluster(8) else cl <- NULL

#Define replicates and groups
replicates <- as.character(metadatah$Groups)
groups <- list(NDE = rep(1,nrow(metadatah)), DE = rep(1:2, each = nrow(metadatah)/2))

#Add annotations
cname <- df[,1]
df = as.matrix(df)

#Create the Countdata object
CD <- new("countData", data = df, replicates = replicates, groups = groups)
dim(df)
dim(CD)

#estimate library sizes for countData object
libsizes(CD) <- getLibsizes(CD)

#Add annotation to the object 
CD@annotation <- as.data.frame(row.names(df))
CD@annotation

#Estimate the prior distributions
CDP.NBML <- getPriors.NB(CD, samplesize = 500, estimation = "QL", cl = cl)

#We then acquire posterior likelihoods, estimating the proportions of differentially expressed counts
CDPost.NBML <- getLikelihoods(CDP.NBML, pET = 'BIC', cl = cl, bootStraps = 3, verbose = F)

#We can ask for the top candidates for differential expression using the topCounts function
dataset2_de = topCounts(CDPost.NBML, group=2, number=Inf)
head(dataset2_de)

# Select based on FDR
dataset2_de <- subset(as.data.frame(dataset2_de), FDR.DE < 0.05)
nrow(dataset2_de)

#Save to file
write.table(dataset2_de, paste0("bayseq",j,".csv"),
            row.names = F,sep=",",col.names=FALSE)
df1=read.csv(paste0("bayseq",j,".csv"),header = F)
source('putinalist.R')
putinalist(df1)



