# Set the working directory
home= "/Volumes/Seagate\ Backup\ Plus\ Drive/PhD/BENIN_git/"
setwd(home)
### libraries for execution of BENIN
library(tseries)
library(xts)
library(boot)
library(caTools)
library(RGeode)
library(dplyr)
library(doSNOW)
library(parallel)
library(foreach)
library(iterators)
library(doSNOW)
library(glmnet)
library(parallel)
library(boot)
library(infotheo)
library(minet)
library(tseries)
library(zoo)
library(xts)
library(boot)
library(caTools)
library(RGeode)
library(crayon)
library(dplyr)
library(data.table)
library(reshape2)
library(parallel)
library(doRNG)
library(data.table)

### 
source("R/utile.R")
source("R/benin.R")
source("src/bootpracs.q")
source("src/bootfuns.q")

inferencenetworkonDREAMNetwork<-function(file,network,lambda,alphaenet,sizenetwork,beta,nbfolds,lambdamin,lambdamax)
{
  #print (lambdamin)
  #rint(lambdamax)
  namefile_exprdata= paste("../data/",
                           "dream4_challenge_data/training_data/",
                           "DREAM4_InSilico_Size",sizenetwork,"/insilico_size",sizenetwork,
                           "_",network,"/",
                           "insilico_size",sizenetwork,"_",network,"_timeseries.tsv",sep="")
  tsdata=read.table(file =namefile_exprdata , sep = '\t', header = TRUE)
  X=as.matrix(tsdata[,-1])
  row.names(X)=tsdata[,1]
  namefile_ref_network=paste("../data/dream4_challenge_data/gold_standard/",
                             "DREAM4_Challenge2_GoldStandards/Size_",sizenetwork,"/",
                             "DREAM4_GoldStandard_InSilico_Size",sizenetwork,"_",network,".tsv",sep="")
  
  gold_standard_network=read.table(file=namefile_ref_network , sep = '\t', header = FALSE)
  #listbootstrap<-c(50,100,500,1000,1500,2000,5000,10000)
  listbootstrap<-seq(0,10000,100)
  listbootstrap<-listbootstrap[-1]
  #listlmean<-seq(5,100,10)
  test=lapply(listbootstrap,FUN= inferencenetwork,tsdata=tsdata,gold_standard_network=gold_standard_network,network=network,lambda=lambda,alphaenet=alphaenet,
         sizenetwork=sizenetwork,beta=beta,nbfolds=10,lambdamin=lambdamin,lambdamax=lambdamax,file=file,lmean=10)
}
inferencenetwork<- function(nbBoobstrap,lmean,file=1,tsdata,gold_standard_network,network,lambda,alphaenet,sizenetwork,beta,nbfolds,lambdamin,lambdamax)
{
  X=as.matrix(tsdata[,-1])
  row.names(X)=tsdata[,1]
  tostringalphaenet=toString(alphaenet)
  tostringalphaenet=gsub(".","",tostringalphaenet,fixed=T)
  #tostringalphaenet
  listtf=as.matrix(unique(gold_standard_network[gold_standard_network$V3==1,]$V1))
  nr=nrow(gold_standard_network)
  nc=ncol(gold_standard_network)
  list_genes=union(gold_standard_network$V1,gold_standard_network$V2)
  numbergenes=length(list_genes)
  
  
  namefilematpval=paste("../data","/dream4_challenge_data/training_data/DREAM4_InSilico_Size",sizenetwork,
                        "/insilico_size",sizenetwork,"_",network,"/lambda_",lambda,"/mat_pval_size",sizenetwork,
                        "_",network,"_lambda_",lambda,"_",file,".tsv",sep = "")
  
  matpval=as.matrix(read.table(file=namefilematpval , sep = '\t', header = TRUE,
                               row.names = 1))
  
  
  
  
  
  
  
  ##create probabilities matrix
  matweight=matrix(0,nr = numbergenes, nc = numbergenes)
  rownames(matweight)<-list_genes
  colnames(matweight)<-list_genes
  for (g1 in list_genes)
  {
    for (g2 in list_genes)
    {
      pval=matpval[g1,g2]
      # print(pval)
      # prob=ifelse(g1==g2,0,(intergral(beta=beta,pval=pval)))
      # matweight[g1,g2]=ifelse(prob>=0.8,0.0001,1)
      matweight[g1,g2]=ifelse(g1==g2,0,1/(integral(beta=beta,pval=pval,lambdamin=lambdamin,lambdamax=lambdamax)))
    }
  }
  #mattf=as.matrix(gold_standard_network$V1[which(gold_standard_network$V3==1)])
  #mattf=apply(mattf,1,function(x){gsub("G","",x,ignore.case = T)})
  #listtf=unique(strtoi(mattf, base = 0L))
  #edgeList=bootstrapbenin(X,nbBoobstrap=nbBoobstrap,matweightpk=matweight,sizenetwork = sizenetwork,indxtf = listtf,nbfolds = nbfolds,alphaenet = alphaenet,lmean=lmean)
  no_cores <- detectCores() - 1
  cl<-makeCluster(no_cores)
  clusterEvalQ(cl, .libPaths("./library")) 
  registerDoSNOW(cl)
  load("random_state_seed1001.RData")
  clusterSetRNGStream(cl)
  #print("iciciiiii")
  edgeList<-bootstrapbenin(X,nbBoobstrap=nbBoobstrap,matweightpk=matweight,sizenetwork = sizenetwork,listtf  = listtf,nbfolds = nbfolds,alphaenet = alphaenet,lmean=lmean,
                           allgenes = list_genes,parallel = TRUE)
  stopCluster(cl)
  global_res<-list2df(edgeList)
  ### Treating the result
  # listgene= names(testboot)
  # global_res=data.frame()
  # for (gene in listgene)
  # {
  #   subres=unlist(testboot[[gene]])
  #   listTF=names(subres)
  #   lenTF=length(subres)
  #   listW=unname(subres)
  #   global_res=rbind(global_res,data.frame(TF=listTF,TG=rep(gene,lenTF),W=listW))
  # }
  test=global_res
  #print(global_res)
  #colnames(global_res)=c("TF","TG","W")
  #global_res[global_res$Regulator==global_res$TargetGene,"W"]=0
  global_res[global_res[,1]==global_res[,2],3]=0
  global_res=global_res[order(global_res$W,decreasing = TRUE),]
  global_res<-global_res%>%subset(TF!=TG)
  #namefileglobalres<-paste("../res_in_sillico",
  #                         "/parameter_checking/var_nbbootstrap/","global_res_",sizenetwork,"_",network,'_',file,"_benin_alphaenet",tostringalphaenet
   #                        ,'_nbbootstrap_',nbBoobstrap,'_lmean_',lmean,'_',sizenetwork,"_",
   #                        network,".tsv",sep = "")
  namefileglobalres<-paste("../res_in_sillico",
                           "/parameter_checking/parameter_checking/",sizenetwork,"_",network,"/var_nbbootstrap/global_res",sizenetwork,"_",network,'_',file,"_benin_alphaenet",tostringalphaenet
                           ,'_nbbootstrap_',nbBoobstrap,'_lmean_',lmean,'_',sizenetwork,"_",
                           network,".tsv",sep = "")
  savedata(data=global_res,file=namefileglobalres,sep="\t",colnames=FALSE, rownames=FALSE,quote=FALSE)
  res_for_eval=data.frame()
  res_for_eval=global_res
  res_for_eval[(res_for_eval$W>=0.5),3]=1
  res_for_eval[(res_for_eval$W<0.5),3]=0
  selected_res_for_eval=as.data.frame(res_for_eval[(res_for_eval$W==1),c(1,2,3)])
  #selected_res_for_eval<-as.data.frame(res_for_eval)
  # 
  # namefilerecNetwork=paste("../res_in_sillico",
  #                          "/insilico_size",sizenetwork,"_",network,"/lambda_",lambda,"/lambda_",lambda,"_alpha_",tostringalphaenet,"/global_res_for_eval_",file,"_",sizenetwork,"_",
  #                          network,".tsv",sep = "")
 
  #namefilerecNetwork<-paste("../res_in_sillico",
 #                         "/parameter_checking/",sizenetwork,"_",network,"/var_nbbootstrap/","global_res_for_eval_",sizenetwork,"_",network,'_',file,"_benin_alphaenet",tostringalphaenet
  #                        ,'_nbbootstrap_',nbBoobstrap,'_lmean_',lmean,'_',sizenetwork,"_",
  #                       network,".tsv",sep = "")
  namefilerecNetwork<-paste("../res_in_sillico",
                           "/parameter_checking/parameter_checking/",sizenetwork,"_",network,"/var_nbbootstrap/","global_res_for_eval_",sizenetwork,"_",network,'_',file,"_benin_alphaenet",tostringalphaenet
                          ,'_nbbootstrap_',nbBoobstrap,'_lmean_',lmean,'_',sizenetwork,"_",
                         network,".tsv",sep = "")
  # namefilerecNetwork=paste("/Users/stephanie_admin/Documents/GRN_inference/",
  #                          "res_in_sillico/insilico_size10_5/lambda_20_alpha_03/global_res_for_eval_10_10_5.tsv",
  #                          sep="")
  #### Network size 100
  # namefilerecNetwork=paste("/Users/stephanie_admin/Documents/GRN_inference/",
  #                          "res_in_sillico/insilico_size100_1/lambda_20_alpha_03/global_res_for_eval_100_1.tsv",
  #                          sep="")
  savedata(data=selected_res_for_eval,file=namefilerecNetwork, colnames=FALSE, rownames=FALSE,sep="\t",quote=FALSE)
}
# listalphaenet<-c(0.1,0.3,0.5,0.7,0.9)
# listnbRun<- c(500,1000,1500,2000,3000,10000)
# listlenblock<-c(5,10,20,30,40,50,60,70)
# fixednbSample<-function (x)
# {
#     listalphaenet<-c(0.1,0.3,0.5,0.7,0.9)
#     listlenblock<-c(5,10,30,50,70)
#     combinaison<-unlist(lapply(listalphaenet, function(X) {lapply(llistlenblock, 
#                                                               function(Y) { c(X, Y) }) }), recursive=FALSE)
# }

#lapply(1:5,inferencenetworkonDREAMNetwork,network=1,lambda=20,alphaenet=0.3,sizenetwork=100,beta=0.5,lambdamax=1000,lambdamin=1)

applyallnetwork<-function(network)
{
  lapply(1:10,inferencenetworkonDREAMNetwork,network=network,lambda=20,alphaenet=0.3,sizenetwork=100,beta=0.5,lambdamax=1000,lambdamin=1)
}
