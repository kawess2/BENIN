#!/usr/bin/env Rscript
setwd("src")
#source("cluster_network_inference.R")
source("human_hela_network_inference.R")
#source("sacc_cellcycle_network_inference.R")
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>=1) {
  # default output file
 # clusterfilenamefile<- paste0("../data/data_human/time_series/result_clustering/finalclustering/",args[1])
 # print (clusterfilenamefile)
  beta<-as.numeric(args[2])
  nbfolds<-as.numeric(args[3])
  lambda<-as.numeric(args[4])
  sizenetwork<-as.numeric(args[5])
  nbBootstrap<-as.numeric(args[6])
  lambdamin<-as.numeric(args[7])
  lambdamax<-as.numeric(args[8])
  lmean<-as.numeric(args[9])
  alphaenet<-as.numeric(args[10])
  exponent<-as.numeric(args[11])
 priordata<-args[12]
#for saccharomyces
methodfuncsim<-args[13]
print(priordata)
print(paste0("the number of args is ",length(args)))
  clusterID<-unlist(strsplit(args[1], split=".",fixed = T))[1]
  testsentence<-paste0("the parameters are beta: ",beta," nbfolds: ",nbfolds," lambda: ",lambda," sizenetwork: ",sizenetwork
                       ," nbBoobstrap: ",nbBootstrap, " lambdamin: ",lambdamin, " lambdamax: ",lambdamax," lmean: ",lmean," alphaenet: ", 
                       alphaenet, " exponent: ", exponent," with prior data", priordata, "the measure for functional prior to compute go sim is ", methodfuncsim)
  print(testsentence)
  #networkinferencecluster(clusterfile=clusterfilenamefile,alphaenet=alphaenet,beta=beta,nbfolds=nbfolds,lmean=lmean,
   #                      sizenetwork=0,nbBoobstrap=nbBoobstrap,lambdamin=lambdamin,lambdamax=lambdamax,exponent=exponent,clusterID=clusterID,methodfuncsim=methodfuncsim)

networkinferencepercluster(priordata=priordata,clusterid=clusterID,lambdamin=lambdamin,lambdamax=lambdamax,beta=beta,exponent=exponent,
                                   alphaenet=alphaenet,nbBootstrap=nbBootstrap,ont="BP",methodfuncsim=methodfuncsim,sizenetwork=sizenetwork)
}



