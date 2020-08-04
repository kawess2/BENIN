home= "/Volumes/Seagate\ Backup\ Plus\ Drive/PhD/BENIN_git/"
setwd(home)


### Libraries for execution of BENIN
library(glmnet)
library(tseries)
library(xts)
library(boot)
library(caTools)
library(RGeode)
library(dplyr)
library('doSNOW')
library('parallel')

#library(data.table,lib.loc="./library")
 library(data.table)
library(Rcpp,lib.loc='library/')
library(AnnotationDbi)
library(org.Sc.sgd.db)
library(parallel)
 
###################
source("dreamfourcode/benin_hyp_param.R")
source("R/utile.R")
source("R/benin.R")
source("R/bootpracs.q")
source("R/bootfuns.q")


inferencenetworkonDREAMNetwork<-function(nbBoobstrap,exponent,file,network,lambda,alphaenet,beta,lambdamin,lambdamax,namefile_exprdata,namefile_ref_network,
                                namefilematpval,lmean=10,nbfolds=10,sizenetwork = 100,potentialTF="all")
{
 
tsdata=read.table(file =namefile_exprdata , sep = '\t', header = TRUE,stringsAsFactors = F)
tsdata=na.omit(tsdata) 
 X=as.matrix(tsdata[,-1])
  row.names(X)=tsdata[,1]
  gold_standard_network=data.table::fread(input=namefile_ref_network , sep = '\t', header = FALSE,stringsAsFactors=F,data.table = F)
  #listbootstrap<-c(50,100,500,1000,1500,2000,5000,10000)
  listbootstrap<-seq(50,10000,100)
 # listlmean<-seq(0,100,10)
#listlmean<-listlmean[-c(1,9)]
listlmean<-lmean
#listalphaenet<-c( seq(0,0.6,0.1),seq(0.8,0.9,0.1))
#listalphaenet<-0.9
#listexponent= seq(0.1,1.5,0.1)
#listexponent<-1
#  listbootstrap<-5000
 # lapply(listbootstrap,FUN= inferencenetwork,alphaenet=listalphaenet,lmean = listlmean,network=network,exponent=listexponent,tsdata=tsdata,gold_standard_network=gold_standard_network,
  #        namefilematpval=namefilematpval, beta=beta,nbfolds=nbfolds,lambdamin=lambdamin,lambdamax=lambdamax,
#	sizenetwork = sizenetwork,lambda = lambda,file=file,potentialTF="all")


  inferencenetwork(lmean=lmean,nbBoobstrap=nbBoobstrap,tsdata=tsdata,gold_standard_network=gold_standard_network,
                   namefilematpval=namefilematpval,exponent=exponent,
                   beta=beta,nbfolds=nbfolds,lambdamin=lambdamin,lambdamax=lambdamax,alphaenet=alphaenet,sizenetwork = sizenetwork,lambda = lambda,file=file,potentialTF=potentialTF,network=network)
   
}
inferencenetwork<- function(nbBoobstrap,alphaenet,exponent,lmean,file,tsdata,gold_standard_network,namefilematpval,lambda,sizenetwork,beta,nbfolds,lambdamin,
                            lambdamax,network,potentialTF="known")
{
  tostringalphaenet=toString(alphaenet)
  tostringalphaenet=gsub(".","",tostringalphaenet,fixed=T)
  tostringexponent=toString(exponent)
  tostringexponent=gsub(".","",tostringexponent,fixed=T)
#namefileglobalres=paste0("../res_in_sillico/insilico_size",sizenetwork,"_",network,"/lambda_",lambda,"/lambda_",lambda,"_alpha_",
 #             tostringalphaenet,"/globares_allgenes_as_tf_nbbootstrap_",nbBoobstrap,"_alphaenet_",tostringalphaenet,"_lmean_",lmean,"_exponent_",tostringexponent,'_',file,"_",sizenetwork,"_",network,".tsv")
#namefilerecNetwork<-paste0("../res_in_sillico/insilico_size",sizenetwork,"_",network,"/lambda_",lambda,"/lambda_",lambda,"_alpha_",
 #                          tostringalphaenet,"/globares_for_eval_allgenes_as_tf_lmean_",lmean,"_exponent_",tostringexponent,"_",file,"_",sizenetwork,"_",network,".tsv")
#namefileglobalres=paste0("../data/hyp_parameters_tuning/",sizenetwork,"_",network,"/res_hyp_param_tuning/lambda_",lambda,"/lambda_",lambda,"_alpha_",
#              tostringalphaenet,"/globares_test_nbbootstrap_",nbBoobstrap,"_alphaenet_",tostringalphaenet,"_lmean_",lmean,"_exponent_",tostringexponent,'_',file,"_",sizenetwork,"_",network,".tsv")
#namefilerecNetwork<-paste0("../data/hyp_parameters_tuning/",sizenetwork,"_",network,"/res_hyp_param_tuning/lambda_",lambda,"/lambda_",lambda,"_alpha_",
#                           tostringalphaenet,"/globares_for_eval_test_lmean_",lmean,"_exponent_",tostringexponent,"_",file,"_",sizenetwork,"_",network,".tsv") 



 X=as.matrix(tsdata[,-1])
  row.names(X)=tsdata[,1]
  mattf=as.matrix(gold_standard_network$V1[which(gold_standard_network$V3==1)])
  listtf<-unique(mattf)
  print(paste0("the network is ",network))

  nr=nrow(gold_standard_network)
  nc=ncol(gold_standard_network)
  list_genes<-union(gold_standard_network$V1,gold_standard_network$V2)
  numbergenes<-length(list_genes)
  if (potentialTF=="known")
  {
   	print("here we selected to work with known genes")
    mattf=as.matrix(gold_standard_network$V1[which(gold_standard_network$V3==1)])
    listtf<-unique(mattf)
    #namefileglobalres=paste0("../res_in_sillico/insilico_size",sizenetwork,"_",network,"/globares_new_nbbootstrap_",nbBoobstrap,"_alphaenet_",tostringalphaenet,"_KOzcorenew_exp$
  namefileglobalres=paste0("../data/hyp_parameters_tuning/",sizenetwork,"_",network,"/res_hyp_param_tuning/lambda_",lambda,"/lambda_",lambda,"_alpha_",
              tostringalphaenet,"/grid_search/globares_test_nbbootstrap_grid_search_",nbBoobstrap,"_alphaenet_",tostringalphaenet,"_lmean_",lmean,"_exponent_",tostringexponent,'_',file,"_",sizenetwork,"_",network,".tsv")
print(namefileglobalres)
}
  else
  {
   print("Here we consider all the genes as TF")
	 listtf<- list_genes
   # namefileglobalres=paste0("../res_in_sillico/insilico_size",sizenetwork,"_",network,"/globares_new_nbbootstrap_allgeneastf_",nbBoobstrap,"_alphaenet_",tostringalphaenet,"_KO$
	namefileglobalres=paste0("../res_in_sillico/tempres_10/",sizenetwork,"_",network,"/globares_test_nbbootstrap_",nbBoobstrap,"_alphaenet_",tostringalphaenet,"_lmean_",lmean,"_exponent_",tostringexponent,'_',file,"_allgenes_",sizenetwork,"_",network,".tsv")
print(namefileglobalres) 

 }
 print(paste("the number of transcription factor", length(listtf),sep = " "))
  if (file.exists(namefilematpval))
  {

    matpval=data.table::fread(input=namefilematpval , sep = '\t', header = TRUE,
                                                       data.table = F,stringsAsFactors=F)
    row.names(matpval)<-matpval$V1
    matpval<-matpval[,-1]
    matpval<-data.matrix(matpval)
  } 
  else
  {
    print("we can\'t find the matpval file")
    matpval<-gold_standard_network
    matpval$V3=as.double( matpval$V3)
    nbtrlink<-dim(gold_standard_network[gold_standard_network$V3==1,])[1]
    nbfllink<-dim(gold_standard_network[gold_standard_network$V3==0,])[1]
    matpval[matpval$V3==1,3]<-rexptr(nbtrlink,lambda,c(0,1))
    matpval[matpval$V3==0,3]<-runif(nbfllink, min = 0, max = 1)
    
    matpval<-reshape2::dcast(matpval,V2~V1,fill=0)
    rownames(matpval)<-matpval[,1]
    matpval<-as.matrix(matpval[,-1])
    data.table::fwrite(as.data.frame(matpval),file=namefilematpval,sep="\t",row.names=T,col.names=T)
  }

  matweight=matrix(0,nr = numbergenes, nc = numbergenes)
  rownames(matweight)<-list_genes
  colnames(matweight)<-list_genes
    for (g1 in list_genes)
   {
      for (g2 in list_genes)
      {
       pval=matpval[g1,g2]
       matweight[g1,g2]=ifelse(g1==g2,0,1/(integral(beta=beta,pval=pval,lambdamin=lambdamin,lambdamax=lambdamax))**exponent)
    }
   }

  execTimeBenin<-system.time(edgeList<-applybootstrapbenin(X,nbBoobstrap=nbBoobstrap,matweightpk=matweight,sizenetwork = sizenetwork,
                                                           listtf = listtf,nbfolds = nbfolds,alphaenet = alphaenet,
lmean=lmean,allgenes = list_genes,normalize=T))[[3]]
  print(paste("Benin took:", execTimeBenin," to proceed a network of size ", sizenetwork,sep = " "))

  global_res<-list2df(edgeList)
 
  test=global_res
  global_res[global_res[,1]==global_res[,2],3]=0
  global_res=subset(global_res,TF!=TG)
  global_res=global_res[order(global_res$W,decreasing = TRUE),]
print(namefileglobalres)
  savedata(data=global_res,file=namefileglobalres,sep="\t", colnames=FALSE, rownames=FALSE,quote=FALSE)

  res_for_eval=global_res
  res_for_eval[(res_for_eval$W>=0.5),3]=1
  res_for_eval[(res_for_eval$W<0.5),3]=0
  selected_res_for_eval=as.data.frame(res_for_eval[(res_for_eval$W==1),c(1,2,3)])
  
  }


applynetworkdream<-function(network,indxfile,lmean=10,beta=0.5,nbfolds=10,nbfile=11,sizenetwork=100,nbBoobstrap=1000,lambdamin=1,lambdamax=1000,lambda=20,exponent=1,potentialTF="all",beginindxbootstrap,endindxbootstrap)
{
if (network==4)
{
	alphaenet=0.9
	#alphaenet=1
}
else{
if(network==3)
{
alphaenet=0.9
	#alphaenet=1
}
else{
if (network==2)
{
alphaenet= 0.9
#alphaenet=1
}
else{
if (network == 1)
{
alphaenet= 0.9
#alphaenet=1
}
else{
alphaenet= 0.9
}
}
}
}
lapply(indxfile:nbfile,FUN=applyfile,lmean=lmean,alphaenet= alphaenet, potentialTF=potentialTF,network=network,beta=beta,nbfolds=nbfolds,sizenetwork=sizenetwork,nbBoobstrap=nbBoobstrap,lambdamin=lambdamin,lambdamax=lambdamax,lambda=lambda,exponent=exponent,beginindxbootstrap=beginindxbootstrap,endindxbootstrap=endindxbootstrap)
}
applyfile<-function(file,network,beta,nbfolds,lambda,sizenetwork,nbBoobstrap,alphaenet,lambdamin,lambdamax,lmean,exponent,potentialTF="all",beginindxbootstrap,endindxbootstrap)
{

namefile_exprdata= paste("../data/hyp_parameters_tuning/",sizenetwork,"_",network,
	"/insilico_size",sizenetwork,"_",network,"_dream4_timeseries.tsv",sep="")

namefilematpval= paste("../data/hyp_parameters_tuning/",sizenetwork,"_",network,
        "/lambda_",lambda,"/mat_pval_size",sizenetwork,
                     "_",network,"_lambda_",lambda,"_",file,".tsv",sep = "")

namefile_ref_network= paste("../data/hyp_parameters_tuning/",sizenetwork,"_",network,
	"/insilico_size",sizenetwork,"_",network,"_goldstandard.tsv",sep="")
#  namefile_exprdata= paste("../data/",
 #                        "dream4_challenge_data/training_data/",
  #                       "DREAM4_InSilico_Size",sizenetwork,"/insilico_size",sizenetwork,
   #                      "_",network,"/",
    #                     "insilico_size",sizenetwork,"_",network,"_timeseries.tsv",sep="")
#namefile_ref_network=paste("../data/dream4_challenge_data/gold_standard/",
 #                          "DREAM4_Challenge2_GoldStandards/Size_",sizenetwork,"/",
  #                         "DREAM4_GoldStandard_InSilico_Size",sizenetwork,"_",network,".tsv",sep="")
#namefilematpval=paste("../data","/dream4_challenge_data/training_data/DREAM4_InSilico_Size",sizenetwork,
 #                     "/insilico_size",sizenetwork,"_",network,"/lambda_",lambda,"/mat_pval_size",sizenetwork,
  #                    "_",network,"_lambda_",lambda,"_",file,".tsv",sep = "")
#namefile_exprdata=paste0("../data/simulated_yeast_data/network_",sizenetwork,"_",network,"/Yeast_",sizenetwork,"-",network,"_dream4_timeseries.tsv")
#namefile_ref_network=paste0("../data/simulated_yeast_data/goldstandard_networks/Yeast_",sizenetwork,"-",network,"_goldstandard.tsv")
#namefileglobalres=paste0("../res_in_sillico/insilico_size",sizenetwork,"_",network,"/lambda_",lambda,"/lambda_",lambda,"_alpha_",
#              tostringalphaenet,"/globares_new_nbbootstrap_",nbBoobstrap,"_alphaenet_",tostringalphaenet,"_",file,"_",sizenetwork,"_",network,".tsv")
#namefilerecNetwork<-paste0("../res_in_sillico/insilico_size",sizenetwork,"_",network,"/lambda_",lambda,"/lambda_",lambda,"_alpha_",
 #                          tostringalphaenet,"/globares_for_eval_new_",file,"_",sizenetwork,"_",network,".tsv")
#namefilematpval=paste0("../data/simulated_yeast_data/network_",sizenetwork,"_",file,"/lambda_",lambda,"/mat_pval_size",sizenetwork,"_",network,
#"_lambda_",lambda,"_",file,".tsv")
#print(namefile_exprdata)x
#inferencenetworkonDREAMNetwork(file=file,lmean=lmean,network=network,lambda=lambda,alphaenet=alphaenet,sizenetwork=sizenetwork,beta=beta,nbfolds=nbfolds,lambdamin=lambdamin,lambdamax=lambdamax,
 #                              namefile_exprdata=namefile_exprdata,namefile_ref_network=namefile_ref_network,
#                               namefilematpval=namefilematpval,nbBoobstrap=nbBoobstrap,potentialTF=potentialTF,exponent=exponent)
gs <- list( nbBoobstrap= seq(beginindxbootstrap,endindxbootstrap,100),
           exponent= seq(0.1,2,0.1)) %>% 
  purrr::cross_df()
pmap(gs, .f=inferencenetworkonDREAMNetwork,file=file,network=network,lambda=lambda,alphaenet=alphaenet,beta=beta,lambdamin=lambdamin,lambdamax=lambdamax,namefile_exprdata=namefile_exprdata,namefile_ref_network=namefile_ref_network,
                                namefilematpval=namefilematpval,lmean=lmean,nbfolds=nbfolds,sizenetwork=sizenetwork,potentialTF=potentialTF)

}
