##
#This code computes for differentially expressed genes using limma
##
#Import package limma
library("limma")
library("edgeR")

#First create a DGEList object using the edgeR package:
dge <- DGEList(counts=df)

#Create a design matrix
design <- cbind("1"=1,"1vs2"=rep(c(1,2), each = nrow(metadatah)/2))

#apply scale normalization to counts, TMM normalization method perform well from comparative studies.
dge <- calcNormFactors(dge)

#In the limma-trend approach, the counts are converted to logCPM values using edgeR’s cpm function:
#The prior count is used here to damp down the variances of logarithms of low counts.
logCPM <- cpm(dge, log=FALSE, prior.count=2)

# The logCPM values can then be used in any standard limma pipeline, 
# using the trend=TRUE argument when running eBayes or treat . 
fit <- lmFit(logCPM, design)
fit <- eBayes(fit)
res_limma_trend_eBayes <- topTable(fit, coef=ncol(design),number=Inf)
res_limma_trend_eBayes <- subset(res_limma_trend_eBayes, adj.P.Val < 0.05)
write.table(as.data.frame(rownames(res_limma_trend_eBayes)),
            paste0("limma-trend_eBayes.csv",j,".csv"),
          row.names = F,sep=",",col.names=FALSE)
df1=read.csv(paste0("limma-trend_eBayes.csv",j,".csv"),header = F)
source('putinalist.R')
putinalist(df1)

# Or, to give more weight to fold-changes in the gene ranking, one might use:
# The logCPM values can then be used in any standard limma pipeline,
fit <- lmFit(logCPM, design)
fit <- treat(fit, lfc = log2(1.2))
res_limma_trend_treat <- topTreat(fit, coef = ncol(design),number=Inf)
res_limma_trend_treat <- subset(res_limma_trend_treat, adj.P.Val < 0.05)
write.table(as.data.frame(rownames(res_limma_trend_treat)),
          paste0("limma-trend_treat.csv",j,".csv"),
          row.names = F,sep=",",col.names=FALSE)
df1=read.csv(paste0("limma-trend_treat.csv",j,".csv"),header = F)
source('putinalist.R')
putinalist(df1)

# The voom transformation is applied to the normalized and filtered DGEList object
#  (When the library sizes are quite variable between samples, 
#    the voom approach is theoretically  more powerful than limma-trend)                              
# The voom transformation uses the experiment design matrix, and produces an EList object.
v <- voom(dge, design)

#It is also possible to give a matrix of counts directly to voom without TMM normalization, by
v2 <- voom(df, design, plot=TRUE)

#If the data are very noisy, one can apply the between-array normalization methods.
v3 <- voom(df, design, plot=TRUE, normalize="quantile")

# After this, the usual limma pipelines for differential expression is be applied.
fit <- lmFit(v, design)
fit <- eBayes(fit)
res_limma_voom_eBayes <- topTable(fit, coef=ncol(design),number=Inf)
res_limma_voom_eBayes <- subset(res_limma_voom_eBayes, adj.P.Val < 0.05)
write.table(as.data.frame(rownames(res_limma_voom_eBayes)),
          paste0("limma_voom_eBayes.csv",j,".csv"),
          row.names = F,sep=",",col.names=FALSE)
df1=read.csv(paste0("limma_voom_eBayes.csv",j,".csv"),header = F)
source('putinalist.R')
putinalist(df1)

#Or, to give more weight to fold-changes in the ranking, one could use say:
fit <- lmFit(v, design)
fit <- treat(fit)
res_limma_voom_treat <- topTreat(fit, coef=ncol(design),number=Inf)
head(res_limma_voom_treat)

res_limma_voom_treat <- subset(res_limma_voom_treat, adj.P.Val < 0.05)
write.table(as.data.frame(rownames(res_limma_voom_treat)),
          paste0("limma_voom_treat.csv",j,".csv"),
          row.names = F,sep=",",col.names=FALSE)
df1=read.csv(paste0("limma_voom_treat.csv",j,".csv"),header = F)
source('putinalist.R')
putinalist(df1)
