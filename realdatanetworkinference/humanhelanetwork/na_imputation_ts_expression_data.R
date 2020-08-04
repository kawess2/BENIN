#library("DMwR")
#library("conclust")
library(curl,lib.loc='library/')
library(DMwR,lib.loc='library/')
setwd("/nfs/speed-scratch/s_kamgni/GRN_inference_prime/src")
human_expr_file<-"../data/data_human/time_series/dataPlusScores_all5.txt"
all_human_expr_data<-read.table(human_expr_file,header=T, stringsAsFactors =F, fill=T,quote="",sep="\t")
all_human_expr_data<-all_human_expr_data[, colSums(is.na(all_human_expr_data)) != nrow(all_human_expr_data)]
all_human_expr_data<-all_human_expr_data[,c(1,13:126)]
row.names(all_human_expr_data)<-all_human_expr_data$UID
all_human_expr_data<-all_human_expr_data[,-1]
print(all_human_expr_data[,1:2])
final_annot_human_exprdata_unique_id<-knnImputation(all_human_expr_data,k=12)
write.table(final_annot_human_exprdata_unique_id, file="../data/data_human/time_series/imputed_data_human_unique_ens_id.txt",row.names = F,quote=F,sep= "\t")
