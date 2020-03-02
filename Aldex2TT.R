##
#This code computes for differentially expressed genes using ALDEx2
##
#Import ALDEx2
library("ALDEx2")
#set conditions
conditions <- as.character(metadatah$Groups)
class(conditions)
# General -----------------------------------------------------------------
Aldy_General <- aldex(df, conditions, mc.samples=128, test="t", effect=TRUE,
                      include.sample.summary=TRUE, denom="all", verbose=TRUE)

##extract differentially expressed genes 
sum(Aldy_General$we.eBH < 0.1)

#######[1] 1175
Aldy_DE1 <- subset(Aldy_General, Aldy_General$we.eBH < 0.1)

#View(Aldy_DEG)
Aldex_DE1 <- rownames(Aldy_DE1)

write.table(as.data.frame(Aldex_DE1),paste0("Aldex",j,".csv"),row.names = F,sep=",",col.names=FALSE)
df1=read.csv(paste0("Aldex",j,".csv"),header = F)
source('putinalist.R')
putinalist(df1)


