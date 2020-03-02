##
#This code computes for differentially expressed genes using edgeR
##
#Loading the edgeR package
#Load the needed libraries
library("edgeR")

#A grouping factor can be added at the same time:
Groups = rep(c(1,2), each = nrow(metadatah)/2)
dge = DGEList(df, group = Groups)

#apply scale normalization to counts, 
#TMM normalization method perform well from comparative studies.
#TMM normalization is applied to this dataset to account for 
#compositional difference between the libraries.
dge <- calcNormFactors(dge)

#Create a design matrix
data.frame(Sample=colnames(dge),Groups)
design <- model.matrix(~Groups)
rownames(design) <- colnames(dge)

#Next we estimate dispersion, to estimate common dispersion and tagwise dispersions we
dge <- estimateDisp(dge,design, robust=F)

#First fit genewise glms:
fit <- glmFit(dge, design)

##the glmRT approach##
#Conduct likelihood ratio tests for 1 vs 2 and show the top genes:
lrt <- glmLRT(fit,coef=2)
#topTags(lrt)

#The total number of differentially expressed genes at 5% FDR is given by:
summary(decideTests(lrt))

#extract all the enes with their logFC, logCPM, LR, PValue and FDR
res_glmRT<-topTags(lrt, n = Inf)
difexpgenes <- subset(as.data.frame(res_glmRT), FDR < 0.05)
nrow(difexpgenes)
write.table(as.data.frame(rownames(difexpgenes)), paste0("EdgeR_glmRT", j,".csv"),
            row.names = F,sep=",",col.names=FALSE)
df1=read.csv(paste0("EdgeR_glmRT", j,".csv"),header = F)
source('putinalist.R')
putinalist(df1)

##the glm quasi-likelihood F-test approach##
fit <- glmQLFit(dge, design)

#Conduct quasi-likelihood F-test for 1 vs 2 and show the top genes:
qlf <- glmQLFTest(fit)
summary(decideTests(qlf))
res_glmQLF<-topTags(qlf, n = Inf)
difexpgenes <- subset(as.data.frame(res_glmQLF), FDR < 0.05)
nrow(difexpgenes)
write.table(as.data.frame(rownames(difexpgenes)),
            paste0("EdgeR_glmQLF",j,".csv"),
            row.names = F,sep=",",col.names=FALSE)
df1=read.csv(paste0("EdgeR_glmQLF",j,".csv"),header = F)
source('putinalist.R')
putinalist(df1)

## the classic EdgeR approach##
dge <- DGEList(counts=df, group=Groups)
dge <- estimateDisp(dge)
et <- exactTest(dge)
summary(decideTests(et))
EdgeR_EXACT<-topTags(et, n = Inf)
difexpgenes <- subset(as.data.frame(EdgeR_EXACT), FDR < 0.05)
nrow(difexpgenes)
write.table(as.data.frame(rownames(difexpgenes)),
            paste0("EdgeR_EXACT",j,".csv"),
            row.names = F,sep=",",col.names=FALSE)
df1=read.csv(paste0("EdgeR_EXACT",j,".csv"),header = F)
source('putinalist.R')
putinalist(df1)