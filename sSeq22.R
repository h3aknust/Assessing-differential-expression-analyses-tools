##
#This code computes for differentially expressed genes using sSeq
##
#Import package sSeq
library(sSeq)

# Set 1 -----------------------------------------------------
conditions <- as.character(metadatah$Groups)

res1 = nbTestSH(df, conditions, "1", "2")
#head(res1)

res1$FDR <- p.adjust(res1$pval, method = "BH")

sSeq_DEG1 <- subset(res1, res1$FDR < 0.05)
sSeq_DEG1 <- rownames(sSeq_DEG1)

disp1 <- nbTestSH(df, conditions, "1", "2", SHonly=TRUE, plotASD=TRUE)

plotDispersion(disp1, legPos="bottomright")

#write.csv(sSeq_DEG1, file = "sSeq_DEG1.csv", row.names = FALSE)
write.table(as.data.frame(sSeq_DEG1),paste0("sSeq",j,".csv"),row.names = F,sep=",",col.names=FALSE)
df1=read.csv(paste0("sSeq",j,".csv"),header = F)
source('putinalist.R')
putinalist(df1)

