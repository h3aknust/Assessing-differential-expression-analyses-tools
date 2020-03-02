##
#This code computes common genes between a tool and the truthset
##

putinalist <- function(df1)
{
  df0=read.csv("Truthset.csv",header = F)
  #df1=read.csv("deseq2sim_data_20.csv",header = F)
  truue=unlist(df0$V1)
  tool=unlist(df1$V1)
  commonones <- intersect(truue,tool); tp <- length(commonones)
  commontotool <-setdiff(tool,commonones); fp <- length(commontotool)
  commontotruue <-setdiff(truue,commonones); fn <- length(commontotruue)
  pre <- tp/(tp+fp); rec <- tp/(tp+fn); f1sc <- (2*pre *rec)/(pre+rec)
  mydata1 <<- c(mydata1,as.character(pre))
  
  mydata2 <<- c(mydata2,as.character(rec))
  
  mydata3 <<- c(mydata3,as.character(f1sc))
  
  mydata4 <<- c(mydata4,as.character(length(tool)))
}
