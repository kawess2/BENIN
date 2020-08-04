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

inferencenetworkonDREAMNetwork<-function(file,network,lambda,alphaenet,beta,lambdamin,lambdamax,namefile_exprdata,namefile_ref_network,namefileglobalres,
                                         namefilerecNetwork,namefilematpval, nbBoobstrap=1000,lmean=10,nbfolds=10,sizenetwork = 100)
{

  tsdata=read.table(file =namefile_exprdata , sep = '\t', header = TRUE,stringsAsFactors = F,blank.lines.skip = TRUE)
  tsdata=na.omit(tsdata)
  X=as.matrix(tsdata[,-1])
  row.names(X)=tsdata[,1]

  gold_standard_network=data.table::fread(input=namefile_ref_network , sep = '\t', header = FALSE,stringsAsFactors=F,data.table = F)

listexponent= seq(0.1,1.5,0.1)
#listexponent=1.4

lapply(listexponent,inferencenetwork,lmean=lmean,nbBoobstrap=nbBoobstrap,tsdata=tsdata,gold_standard_network=gold_standard_network,
                                     namefilematpval=namefilematpval,namefileglobalres=namefileglobalres,namefilerecNetwork=namefilerecNetwork,network=network,
                                    beta=beta,nbfolds=nbfolds,lambdamin=lambdamin,lambdamax=lambdamax,alphaenet=alphaenet,sizenetwork = sizenetwork,lambda = lambda,file=file)
   
}
inferencenetwork<- function(exponent,lmean,nbBoobstrap,file,tsdata,gold_standard_network,namefilematpval,lambda,alphaenet,sizenetwork,beta,nbfolds,lambdamin,
                            lambdamax,namefileglobalres,namefilerecNetwork,network)
{
  tostringalphaenet=toString(alphaenet)
  tostringalphaenet=gsub(".","",tostringalphaenet,fixed=T)
  tostringexponent=toString(exponent)
  tostringexponent=gsub(".","",tostringexponent,fixed=T)
  namefileglobalres=paste0("../res_in_sillico/insilico_size",sizenetwork,"_",network,"/lambda_",lambda,"/lambda_",lambda,"_alpha_",
                           tostringalphaenet,"/globares_new_nbbootstrap_",nbBoobstrap,"_exponent_",tostringexponent,"_alphaenet_",tostringalphaenet,"_",file,"_both_",sizenetwork,"_",network,".tsv")
  namefilerecNetwork<-paste0("../res_in_sillico/insilico_size",sizenetwork,"_",network,"/lambda_",lambda,"/lambda_",lambda,"_alpha_",
                             tostringalphaenet,"/globares_for_eval_nbbootstrap2_",nbBoobstrap,"_exponent_",tostringexponent,"_alphaenet_",tostringalphaenet,"_",file,"_both_",sizenetwork,"_",network,".tsv")

  print(namefileglobalres)
  X=as.matrix(tsdata[,-1])
  row.names(X)=tsdata[,1]
  mattf=as.matrix(gold_standard_network$V1[which(gold_standard_network$V3==1)])

  listtf<-unique(mattf)
  print(paste("the number of transcription factor", length(listtf),sep = " "))
  nr=nrow(gold_standard_network)
  nc=ncol(gold_standard_network)
  allgenes=union(gold_standard_network$V1,gold_standard_network$V2)
  numbergenes=length(allgenes)
  
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
    data.table::fwrite(matpval,file=namefilematpval,sep="\t")
  }
 

  ##create probabilities matrix
  matweightlocdata=matrix(0,nr = numbergenes, nc = numbergenes)
  rownames(matweightlocdata)<-allgenes
  colnames(matweightlocdata)<-allgenes
   for (g1 in allgenes)
  {
     for (g2 in allgenes)
     {
      pval=matpval[g1,g2]
      matweightlocdata[g1,g2]=ifelse(g1==g2,0,1/(integral(beta=beta,pval=pval,lambdamin=lambdamin,lambdamax=lambdamax))**exponent)
   }
   }
  kodatafile=paste("../data/",
                   "dream4_challenge_data/training_data/",
                   "DREAM4_InSilico_Size",sizenetwork,"/insilico_size",sizenetwork,
                   "_",network,"/",
                   "insilico_size",sizenetwork,"_",network,"_knockouts.tsv",sep="")
  matkodata<-read.table(file =kodatafile , sep = '\t', header = TRUE,stringsAsFactors = F)
  rownames(matkodata)<-colnames(matkodata)
  wtexprdatafile<-paste("../data/",
                        "dream4_challenge_data/training_data/",
                        "DREAM4_InSilico_Size",sizenetwork,"/insilico_size",sizenetwork,
                        "_",network,"/",
                        "insilico_size",sizenetwork,"_",network,"_wildtype.tsv",sep="")
  wtexprdata<-read.table(file =wtexprdatafile , sep = '\t', header = TRUE,stringsAsFactors = F)
  matkodata_all<-rbind(matkodata,wtexprdata)
  vectsd=apply(matkodata_all, 2, sd)
  vectmean=apply(matkodata_all, 2, mean)
  zscore=apply(matkodata,1,FUN=function(x,vectmean,vectsd){(x-vectmean)/vectsd},vectmean,vectsd)
  matpval=apply(zscore,2,FUN=function(x){2*pnorm(as.matrix(-abs(x)))})
  rownames(matpval)<-colnames(matpval)
  matweightKO=matrix(0,nr = numbergenes, nc = numbergenes)
  rownames(matweightKO)<-allgenes
  colnames(matweightKO)<-allgenes
  for (g1 in allgenes)
  {
    for (g2 in allgenes)
    {
      pval=matpval[g1,g2]
      matweightKO[g1,g2]=ifelse(g1==g2,0,1/(integral(beta=beta,pval=pval,lambdamin=lambdamin,lambdamax=lambdamax))**exponent)
    }
  }
  listprior<- list(matweightKO,matweightlocdata)

colnames(gold_standard_network)<-c("TF","TG","connectivity")
dataframefromlodata<-as.data.frame(as.table(matweightlocdata),stringsAsFactors = F)
dataframefromKO<-as.data.frame(as.table(matweightKO),stringsAsFactors = F)

colnames(dataframefromKO)<-c("TG","TF","W_KO")
colnames(dataframefromlodata)<-c("TG","TF","W_locdata")
mydata<-full_join(gold_standard_network,dataframefromlodata,by=c("TF","TG"))
mydata<-full_join(mydata,dataframefromKO,by=c("TF","TG"))
mydata[is.na(mydata$connectivity),"connectivity"]<-0
logitmodel <- glm(connectivity ~ W_locdata + W_KO, data = mydata, family = "binomial")
mydata["combinedweight"] <-predict(logitmodel,type="response")
   matweightpk <- reshape2::dcast(mydata[,c("TF","TG","combinedweight")],TG~TF,fill=0)
    rownames(matweightpk)<-matweightpk[,1]
    matweightpk<-as.matrix(matweightpk[,-1])

matweightpk=1./(matweightpk**exponent)
matweightpk[is.infinite(matweightpk)]<-0
print(summary(logitmodel))
  execTimeBenin<-system.time(edgeList<-applybootstrapbenin(X,nbBoobstrap=nbBoobstrap,matweightpk=matweightpk,sizenetwork = sizenetwork,
                                                           listtf = listtf,nbfolds = nbfolds,alphaenet = alphaenet,lmean=lmean,allgenes = allgenes))[[3]]
  print(paste("Benin took:", execTimeBenin," to proceed a network of size ", sizenetwork,sep = " "))

 
  global_res<-list2df(edgeList)
  global_res$TF<-as.character(global_res$TF)
  global_res$TG<-as.character(global_res$TG)
  test=global_res
  global_res[global_res[,1]==global_res[,2],3]=0
  global_res=subset(global_res,TF!=TG)
  global_res=global_res[order(global_res$W,decreasing = T),]
  
  savedata(data=global_res,file=namefileglobalres,sep="\t", colnames=FALSE, rownames=FALSE,quote=FALSE)

  res_for_eval=global_res
  res_for_eval[(res_for_eval$W>=0.5),3]=1
  res_for_eval[(res_for_eval$W<0.5),3]=0
  selected_res_for_eval=as.data.frame(res_for_eval[(res_for_eval$W==1),c(1,2,3)])

  
  savedata(data=selected_res_for_eval,file=namefilerecNetwork, colnames=FALSE, rownames=FALSE,sep="\t",quote=FALSE)
}


applynetworkdream<-function(network,indxfile,lmean=10,beta=0.5,nbfolds=10,nbfile=11,sizenetwork=100,nbBoobstrap=1000,lambdamin=1,lambdamax=1000,lambda=20)
{
	if (network==4)
	{
        	alphaenet=0.8
       
	}
	else{
		if(network==3)
		{
        		alphaenet=0.7
		}
	else{
	if (network==2)
	{
		alphaenet= 0.8

	}
	else{
	if (network == 1)
	{
	alphaenet= 0.9

	}
	else{
	alphaenet= 0.9
	}
	}
	}
	}




lapply(indxfile:nbfile,FUN=applyfile,lmean=lmean,alphaenet= alphaenet, network=network,beta=beta,nbfolds=nbfolds,sizenetwork=sizenetwork,nbBoobstrap=nbBoobstrap,lambdamin=lambdamin,lambdamax=lambdamax,lambda=lambda)
}






applyfile<-function(file,network,lmean,alphaenet,beta,nbfolds,lambda,sizenetwork,nbBoobstrap,lambdamin,lambdamax)
{
  tostringalphaenet=toString(alphaenet)
  tostringalphaenet=gsub(".","",tostringalphaenet,fixed=T)
  namefile_exprdata= paste("../data/",
                         "dream4_challenge_data/training_data/",
                         "DREAM4_InSilico_Size",sizenetwork,"/insilico_size",sizenetwork,
                         "_",network,"/",
                         "insilico_size",sizenetwork,"_",network,"_timeseries.tsv",sep="")
namefile_ref_network=paste("../data/dream4_challenge_data/gold_standard/",
                           "DREAM4_Challenge2_GoldStandards/Size_",sizenetwork,"/",
                           "DREAM4_GoldStandard_InSilico_Size",sizenetwork,"_",network,".tsv",sep="")
namefilematpval=paste("../data","/dream4_challenge_data/training_data/DREAM4_InSilico_Size",sizenetwork,
                      "/insilico_size",sizenetwork,"_",network,"/lambda_",lambda,"/mat_pval_size",sizenetwork,
                      "_",network,"_lambda_",lambda,"_",file,".tsv",sep = "")#namefile_ref_network=paste0("../data/simulated_yeast_data/goldstandard_networks/Yeast_",sizenetwork,"-",network,"_goldstandard.tsv")
namefileglobalres=paste0("../res_in_sillico/insilico_size",sizenetwork,"_",network,"/lambda_",lambda,"/lambda_",lambda,"_alpha_",
              tostringalphaenet,"/globares_new_nbbootstrap_",nbBoobstrap,"_alphaenet_",tostringalphaenet,"_",file,"_",sizenetwork,"_",network,".tsv")
namefilerecNetwork<-paste0("../res_in_sillico/insilico_size",sizenetwork,"_",network,"/lambda_",lambda,"/lambda_",lambda,"_alpha_",
                           tostringalphaenet,"/globares_for_eval_nbbootstrap2_",nbBoobstrap,"_alphaenet_",tostringalphaenet,"_",file,"_",sizenetwork,"_",network,".tsv")


inferencenetworkonDREAMNetwork(file=file,network=network,lambda=lambda,alphaenet=alphaenet,sizenetwork=sizenetwork,beta=beta,nbfolds=nbfolds,lambdamin=lambdamin,lambdamax=lambdamax,
                               namefile_exprdata=namefile_exprdata,namefile_ref_network=namefile_ref_network,namefileglobalres=namefileglobalres,namefilerecNetwork=namefilerecNetwork,
                               namefilematpval=namefilematpval,lmean=10,nbBoobstrap=nbBoobstrap)
}


