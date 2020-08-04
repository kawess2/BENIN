######## Commented package to load ##################
#library(glmnet)
#library(tseries)
#library(xts)
#library(boot)
#library(caTools)
#library(RGeode)
#library(dplyr)
#library('doSNOW')
#library('parallel')
#library('foreach')
#library('snow')
#library('iterators')
#library('doSNOW')
#library('glmnet')
#library('parallel')
#library('boot')
#library(tseries)
#library(xts)
#library(boot)
#library(caTools)
#library(RGeode)
#library(crayon)
#library(dplyr)
#library(data.table)
#library(reshape2)
#library('parallel')
#library(doRNG)
##################################################################
library(glmnet)
library(boot)
library(dplyr)
library('doSNOW')
library('parallel')
library('foreach')
library('snow')
library('iterators')
library('doSNOW')
library(crayon)
library(dplyr)
library(data.table)
library(reshape2)
library('parallel')
library(doRNG)
library('parallel')

#######
home= "/Volumes/Seagate\ Backup\ Plus\ Drive/PhD/BENIN_git/"
setwd(home)
source('R/Functions_Valid.R')
source("R/bootpracs.q")
source("R/bootfuns.q")
source("R/list2df.R")
source("R/benin2.R")
source("beninapplication/realdatanetworkinference/humanhelanetwork/human_hela_network_inference.R")






networkinferencecluster<-function(clusterfile,alphaenet,beta,nbfolds,sizenetwork=0,nbBoobstrap,lambdamin,lambdamax,lmean,exponent,clusterID){

tostringalphaenet=toString(alphaenet)
  tostringalphaenet=gsub(".","",tostringalphaenet,fixed=T)
  tostringexponent=toString(exponent)
  tostringexponent=gsub(".","",tostringexponent,fixed=T)



  imputed_annot_human_exprdata_unique_id<-read.table(file="../data/data_human/time_series/imputed_ts_data_human_unique_ens_id.txt",
                                                     sep= "\t",stringsAsFactors = F,header=T,row.names=1)
  
  namecluster<-strsplit(clusterfile,"/")[[1]][5]
  imputed_annot_human_exprdata_unique_id=scale(x=imputed_annot_human_exprdata_unique_id,scale = apply(imputed_annot_human_exprdata_unique_id, 2,
                                                                                                      sd, na.rm = TRUE))
listgenesWithExprData <- rownames(imputed_annot_human_exprdata_unique_id)  
listgenes=read.table(file=clusterfile,sep= "\t",stringsAsFactors = F,header=F)

  listgenes<-intersect(listgenes$V1, listgenesWithExprData)
  infoalltfs<-read.csv("../data/data_human/List_TF/DatabaseExtract_v_1.01.txt",header = T,stringsAsFactors = F,sep = "\t")
  infoalltfs<-subset(infoalltfs,Is.TF.=="Yes")
  listtf<-unique(infoalltfs$Ensembl.ID)
  listtfKD<-read.table(file="../data/data_human/geneKnockdown/listKoTFEnsID.txt",sep= "\t",stringsAsFactors = F,header=F)

listtf<-unique(listtfKD$V1) 
listtf<-intersect(listtf,listgenesWithExprData)  
matkddata<-read.table(file="../data/data_human/geneKnockdown/copymatKDpval.txt",sep= "\t",stringsAsFactors = F,header=T,row.names = 1)
  
  listgenes<-intersect(listgenes,unique(rownames(matkddata)))
  allgenes<-union(listgenes,listtf)
  nbtf<-length(listtf)
  numbergenes<-length(allgenes)
  matweight=matrix(0,nr = numbergenes, nc = numbergenes)
  rownames(matweight)<-allgenes
  colnames(matweight)<-allgenes
  matkddata[is.na(matkddata)]<-max(matkddata[!is.na(matkddata)])

for (g1 in allgenes)
  {
    for (g2 in listtf)
    {
	if (g1 %in%unique(rownames(matkddata)) & g2 %in%unique(colnames(matkddata)))
		{
      pval=matkddata[g1,g2]
    	matweight[g1,g2]=ifelse(g1==g2,0,1/(integral(beta=beta,pval=pval,lambdamin=lambdamin,lambdamax=lambdamax))**exponent)
    }
}
  }

matweight[matweight==0]<-max(matweight)
print(matweight)
print("Done computing the matweight matrix")
matexprdata<-as.matrix(t(imputed_annot_human_exprdata_unique_id[allgenes,]))
print(head(matexprdata))
execTimeBenin<-system.time(edgeList<-applybootstrapbenin(matexprdata,nbBoobstrap=nbBoobstrap,matweightpk=matweight,sizenetwork = sizenetwork,
                                                         listtf = listtf,nbfolds = nbfolds,alphaenet = alphaenet,lmean=lmean,allgenes = allgenes))[[3]]
print(paste("Benin took:", execTimeBenin," to proceed a network of size ", sizenetwork,sep = " "))
global_res<-list2df(edgeList)
test=global_res
global_res[global_res[,1]==global_res[,2],3]=0
global_res=subset(global_res,TF!=TG)
global_res=global_res[order(global_res$W,decreasing = TRUE),]

#change the name
#We create a folder to store the results if it doesn't exist already
dir.create("../res_in_sillico/human_cluster_networkinference/",showWarnings = FALSE)
namefileglobalres<-paste0("../res_in_sillico/human_cluster_networkinference/",clusterID,"_global_res_alphaenet_",tostringalphaenet
                         ,'_nbbootstrap_',nbBoobstrap,'_lmean_',lmean, ".tsv")
savedata(data=global_res,file=namefileglobalres,sep="\t",colnames=NA)

}
