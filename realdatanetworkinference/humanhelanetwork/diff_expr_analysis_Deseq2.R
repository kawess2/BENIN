setwd("/Users/stephanie_admin/Documents/GRN_inference/src/")
library(DESeq2)

#Reading count matrix
countfile<-"../data/data_human/final_data_hum_reg_network/Hela_data/KD_data/KD_CTCF/GSE108869_Knockdown_RNA-seq_ensID.txt"
countmat<-read.table(countfile,header = T,stringsAsFactors = F,quote = "",sep = "\t",row.names = 1)
sampleinfile<-"../data/data_human/final_data_hum_reg_network/Hela_data/KD_data/KD_CTCF/info_sample.txt"
infosample<-read.table(sampleinfile,header = T,stringsAsFactors = F,quote = "",sep = "\t",row.names = 1)
infosample$condition<-factor(infosample$condition,levels = unique(infosample$condition))
infosample<-infosample[colnames(countmat),]
ddsMat<-DESeqDataSetFromMatrix(countData = countmat,colData = infosample,design = ~condition)
dds<-DESeq(ddsMat)
#res<-results(dds,contrast = c("condition", "CtrlshRNA","shCTCF"),lfcThreshold = 0.5,alpha=0.05)
res<-lfcShrink(dds,contrast = c("condition", "CtrlshRNA","shCTCF"))