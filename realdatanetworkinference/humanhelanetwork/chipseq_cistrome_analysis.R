library(data.table)
home="Volumes/Seagate\ Backup\ Plus\ Drive/PhD/GRN_inference_prime/src"
setwd(home)
cistromefolder="../data/data_human/final_data_hum_reg_network/Hela_data/chipseqdata/wgEncodeAwgTfbsUniform/cistrome_output"
cistromeinfofileres=paste0(cistromefolder,"/cistroresdatainfo.txt")
cistrome<- read.table(cistromeinfofileres,header = T,stringsAsFactors = F, quote = "",sep="\t")
lisresfile<-cistrome$filerescistrome
finalres<-data.frame()
for (file in lisresfile)
{
  fileres<-paste0(cistromefolder,"/",file)
  res<-fread(fileres,header = T,stringsAsFactors = F, quote = "",sep="\t", skip=4,
             colClasses = c("character","numeric","numeric","character","numeric","character","character"),data.table=FALSE)
  currTF<-subset(cistrome,cistrome$filerescistrome==file,select=c("TF"))$TF
  nrows<-dim(res)[1]
  res$TF<-rep(currTF,nrows)
 
  finalres<-rbind(finalres,res[,c("TF","GeneSymbol","Score")])
}
write.table(finalres,"../data/data_human/final_data_hum_reg_network/Hela_data/chipseqdata/wgEncodeAwgTfbsUniform/TFTGscore_chipseq.txt",
            sep="\t",col.names = T,row.names = F, quote=F)