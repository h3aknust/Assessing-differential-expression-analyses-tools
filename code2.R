#
#
#This code calls all tools used in this analysis and creates csv files for 
#1. selected enes
#2. precision
#3. recall 
#4. fscore
#
#


#create 3 dataframes
mydata0 <-c("Dataset","baySeq.R","deseq","edger","edger","edger","limma","limma","limma","limma",'Absseq.R','AbsseqAFold.R','Aldex2TT.R','ebseq.R','sSeq22.R')
test2 <- list(mydata0)
df00 <- as.data.frame(test2)
colnames(df00) <- df00[1,]
myfscore <- df00
myprecision <- df00
myrecall <- df00
myselected <- df00
#source(putinalist.R)
namelist <- c()
for(i in seq(10,15,5)) {
  mynames <- paste0("simdata_", i)
  namelist <- c(namelist, mynames) 
}
#print(namelist)
for(j in namelist) {
  mydata1 <- c()
  mydata2 <- c()
  mydata3 <- c()
  mydata4 <- c()
  mydata1 <- c(mydata1,j)
  mydata2 <- c(mydata2,j)
  mydata3 <- c(mydata3,j)
  mydata4 <- c(mydata4,j)
  
  #Create names for the simulated data
  source('Simulate.R')
  
  #Read the count data into R dataframe
  df=read.csv(paste0(j,".csv"))
  rownames(df) <- df[,1]
  df <- df[,-1]
  
  ## Upload metadatah
  metadatah=read.csv(paste0(mynames,"_metadata.csv"))
  metadatah=read.csv(paste0(j,"_meta.csv"))
  source('baySeq.R')
  source('DESeq2.R')
  source('edgeR.R')
  source('limma.R')
  source('Absseq.R')
  source('Aldex2TT.R')
  source('ebseq.R')
  source('sSeq22.R')
  
  #Create rows to be used as column names...
  test2 <- list( mydata1)
  df10 <- as.data.frame(test2)
  colnames(df10) <- df10[1,]
  
  #Create table of precision values.
  myprecision <- cbind(myprecision,df10)
  test2 <- list(mydata2)
  df20 <- as.data.frame(test2)
  colnames(df20) <- df20[1,]
  
  #Create table of precision values
  myrecall <- cbind(myrecall,df20)
  test2 <- list(mydata3)
  df30 <- as.data.frame(test2)
  colnames(df30) <- df30[1,]
  
  #Create table of F1-Score values
  myfscore <- cbind(myfscore,df30)
  test2 <- list(mydata4)
  df40 <- as.data.frame(test2)
  colnames(df40) <- df40[1,]
  
  #Create table of Selected genes values
  myselected <- cbind(myselected,df40)
  print("we are done wit round")
  print(j)
  
  #Write perfomance measures to file...

  write.table(as.data.frame(myprecision),"precission.csv",row.names = F,sep=",",col.names=FALSE)

  write.table(as.data.frame(myrecall),"recall.csv",row.names = F,sep=",",col.names=FALSE)

  write.table(as.data.frame(myfscore),"fscore.csv",row.names = F,sep=",",col.names=FALSE)

  write.table(as.data.frame(myselected),"selected.csv",row.names = F,sep=",",col.names=FALSE)
}
