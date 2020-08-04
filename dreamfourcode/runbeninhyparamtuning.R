#!/usr/bin/env Rscript
setwd("src")
source("script_DREAM_Network_inference_mod.R")

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}else if (length(args)>=1) {

indxfile<-args[1]
lmean<-args[2]
beta<-args[3]
nbfolds<-args[4]
nbfile<-args[5]
sizenetwork<-args[6]
nbBoobstrap<-args[7]
lambdamin<-args[8]
lambdamax<-args[9]
lambda<-args[10]
exponent<-args[11]
potentialTF<-args[12]
beginindxbootstrap<-args[13]
endindxbootstrap<-args[14]

print(paste0("The parameters are indxfile: ",indxfile," lmean: ",lmean," beta: ",beta," nbfolds: ",nbfolds," nbfile: ",nbfile," sizenetwork: ",sizenetwork," nbBoobstrap: ",nbBoobstrap,
	" lambdamin: ",lambdamin," lambda: ",lambda," exponent: ", exponent, " potentialTF: ", potentialTF," beginindxbootstrap: ",beginindxbootstrap," endindxbootstrap: ",endindxbootstrap))

lapply(5:5,applynetworkdream,indxfile=indxfile,lmean=lmean,beta=beta,nbfolds=nbfolds,nbfile=nbfile,sizenetwork=sizenetwork,nbBoobstrap=nbBoobstrap,
lambdamin=lambdamin,lambdamax=lambdamax,lambda=lambda,exponent=exponent,potentialTF="know",beginindxbootstrap=beginindxbootstrap,endindxbootstrap=endindxbootstrap)
}
