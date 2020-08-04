setwd("/speed-scratch/s_kamgni/GRN_inference_prime/src")
#suppressMessages(library("universalmotif"))
suppressMessages(library("DMwR"))
suppressMessages(library("conclust"))

save_elmt_cluster<-function(cluster_id,clusters_pred)
{
  namefilecluster<-paste0("../data/data_human/time_series/result_clustering/cluster_mpckm_nbsmpl_500_",cluster_id,".txt")
  genes_in_cluster<-as.data.frame(names(which(clusters_pred==cluster_id)))
  write.table(genes_in_cluster,namefilecluster,sep = "\t",col.names = F,quote = F,row.names = F)
}

#Clustering of TS expression data
nbclusters=35
thres=0.5
nbsubsample=500
# Preparing for clustering. We will use different clustering and combine their idea to get an ensemble clustering
#We are reading the imputed ts data from expression data file
imputed_annot_human_exprdata_unique_id<-read.table(file="../data/data_human/time_series/imputed_ts_data_human_unique_ens_id.txt",
                                                   sep= "\t",stringsAsFactors = F,header=T,row.names=1)


#imputed_annot_human_exprdata_unique_id<-imputed_annot_human_exprdata_unique_id[1:50,]

imputed_annot_human_exprdata_unique_id=scale(x=imputed_annot_human_exprdata_unique_id,scale = apply(imputed_annot_human_exprdata_unique_id, 2, sd, na.rm = TRUE))
listgenes<-rownames(imputed_annot_human_exprdata_unique_id)
#print(head(imputed_annot_human_exprdata_unique_id))
#We load functional annotation data
funcoup<-read.table("../data/data_human/funcoup/FC4.0_H.sapiens_full",header=T,sep = "\t", stringsAsFactors = F,quote="")
subfuncoup<-subset(funcoup,(X2.Gene1%in%listgenes & X3.Gene2%in%listgenes))
mustLink<-subset(subfuncoup,X0.PFC>=thres)
cantLink<-subset(subfuncoup,X0.PFC<thres)
lstmustLink<-mustLink[sample(nrow(mustLink), nbsubsample),c("X2.Gene1","X3.Gene2")]
lstcantLink<-cantLink[sample(nrow(cantLink), nbsubsample),c("X2.Gene1","X3.Gene2")]
lstmustLink[,"Gene1.indx"]<-match(lstmustLink$X2.Gene1, rownames(imputed_annot_human_exprdata_unique_id))
lstmustLink[,"Gene2.indx"]<-match(lstmustLink$X2.Gene2, rownames(imputed_annot_human_exprdata_unique_id))
lstcantLink[,"Gene1.indx"]<-match(lstcantLink$X2.Gene1, rownames(imputed_annot_human_exprdata_unique_id))
lstcantLink[,"Gene2.indx"]<-match(lstcantLink$X3.Gene2, rownames(imputed_annot_human_exprdata_unique_id))

lstmustLink<-lstcantLink[,c("Gene1.indx","Gene2.indx")]
lstcantLink<-lstcantLink[,c("Gene1.indx","Gene2.indx")]
clusters_pred = mpckm(imputed_annot_human_exprdata_unique_id, k=nbclusters, lstmustLink, lstcantLink)
names(clusters_pred )<-rownames(imputed_annot_human_exprdata_unique_id)
listclusters<-unique(clusters_pred)
lapply(listclusters,save_elmt_cluster,clusters_pred=clusters_pred)
