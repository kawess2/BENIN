setwd("/Users/stephanie_admin/Documents/GRN_inference/src/")
library(TFTargetCaller)
library(org.Hs.eg.db)
library(dplyr)
library(reshape2)
mypeakfile<-list.files("../data/data_human/final_data_hum_reg_network/hela_narrow_peak")
mypeakfile<-mypeakfile[-1]
fileinfochipseq<-"../data/data_human/final_data_hum_reg_network/hela_narrow_peak/filenames.csv"
infochipseq<-read.table(fileinfochipseq,stringsAsFactors = F,header=T,quote = "",sep = ",")
subinfochipseq<-subset(infochipseq,infochipseq$filename%in%mypeakfile)
# We get the the gene Position from biomart
genePosition <- getGenePosition("Human", type="protein_coding")
nbtf<-length(subinfochipseq$name)
final_TF_target_score<-as.list(1:nbtf)
names(final_TF_target_score)<-subinfochipseq$name
for (file in mypeakfile)
{
  print(file)
  currentTF<-subset(subinfochipseq,subinfochipseq$filename==file)$name
  if (!length(currentTF)==0)
  {
    print(currentTF)
    filenamepeak<-paste0("../data/data_human/final_data_hum_reg_network/hela_narrow_peak/",file)
    peakPosition <- getPeakPosition(filename=filenamepeak, intensity=7)
    TF_target_score<-TFTargetCaller(peakPosition, genePosition, method="Ouyang", n=5000)
    #print(TF_target_score)
    final_TF_target_score[[currentTF]]<-TF_target_score
  }

  
}

TF_TG_score<-list2df(final_TF_target_score)
colnames(TF_TG_score)<-c("TG","TF","W")
TF_TG_score<-TF_TG_score[,c("TF","TG","W")]
TF_TG_score$TF<-as.character(TF_TG_score$TF)
TF_TG_score$TG<-as.character(TF_TG_score$TG)
mapping_TG_entrez_name<-select(org.Hs.eg.db,unique(TF_TG_score$TG),
                            c("ENTREZID","GENENAME","ENSEMBL","SYMBOL"), "ENSEMBL")
final_TF_TG_score<-inner_join( TF_TG_score,mapping_TG_entrez_name[,c("ENSEMBL","SYMBOL")], by=c("TG"="ENSEMBL"))
final_TF_TG_score<-final_TF_TG_score[,c("TF","SYMBOL","W")]
colnames(final_TF_TG_score)<-c("TF","TG","W")