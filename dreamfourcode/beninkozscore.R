

### Loading libraries
library("BiocGenerics")
library("backports")
library("pkgconfig")

############
source("R/utile.R")
source("script_DREAM_Network_inference_all_prior.R")
source("R/benin.R")
source("src/bootfuns.q")
source("src/bootpracs.q")


##############
BENINKozscore<-function(network=1,alphaenet=0.9,nbBoobstrap=4000,exponent=1.4,sizenetwork=100,lambda=20,lmean=10,alpha=0.5,lambdamin=1,lambdamax=1000,potentialtf="all")
{
  tostringalphaenet=toString(alphaenet)
  tostringalphaenet=gsub(".","",tostringalphaenet,fixed=T)
tostringexponent=toString(exponent)
tostringexponent=gsub(".","",tostringexponent,fixed=T)
  namefile_exprdata= paste("../data/","dream4_challenge_data/training_data/","DREAM4_InSilico_Size",sizenetwork,"/insilico_size",sizenetwork,"_",network,"/","insilico_size",sizenetwork,"_",network,"_timeseries.tsv",sep="")
  print(namefile_exprdata)
  namefile_ref_network=paste("../data/dream4_challenge_data/gold_standard/","DREAM4_Challenge2_GoldStandards/Size_",sizenetwork,"/","DREAM4_GoldStandard_InSilico_Size",sizenetwork,"_",network,".tsv",sep="")
  print(namefile_ref_network)
  namefileglobalres=paste0("../res_in_sillico/insilico_size",sizenetwork,"_",network,"/lambda_",lambda,"/lambda_",lambda,"_alpha_",tostringalphaenet,"/globares_new_nbbootstrap_",nbBoobstrap,"_alphaenet_",tostringalphaenet,"_KO_sumbeta_",sizenetwork,"_",network,".tsv")
  tsdata=read.table(file =namefile_exprdata , sep = '\t', header = TRUE,stringsAsFactors = F)
  tsdata=na.omit(tsdata)
  kodatafile=paste("../data/","dream4_challenge_data/training_data/","DREAM4_InSilico_Size",sizenetwork,"/insilico_size",sizenetwork,"_",network,"/","insilico_size",sizenetwork,"_",network,"_knockouts.tsv",sep="")
  print(kodatafile)
  matkodata<-read.table(file =kodatafile , sep = '\t', header = TRUE,stringsAsFactors = F)
  rownames(matkodata)<-colnames(matkodata)
  wtexprdatafile<-paste("../data/","dream4_challenge_data/training_data/","DREAM4_InSilico_Size",sizenetwork,"/insilico_size",sizenetwork,"_",network,"/",
                        "insilico_size",sizenetwork,"_",network,"_wildtype.tsv",sep="")
  wtexprdata<-read.table(file =wtexprdatafile , sep = '\t', header = TRUE,stringsAsFactors = F)
  print(wtexprdatafile)
  matkodata_all<-rbind(matkodata,wtexprdata)
  vectsd=apply(matkodata_all, 2, sd)
  vectmean=apply(matkodata_all, 2, mean)
  zscore=apply(matkodata,1,FUN=function(x,vectmean,vectsd){(x-vectmean)/vectsd},vectmean,vectsd)
print(zscore)  
gold_standard_network=data.table::fread(input=namefile_ref_network , sep = '\t', header = FALSE,stringsAsFactors=F,data.table = F)
  X=as.matrix(tsdata[,-1])
  row.names(X)=tsdata[,1]
  nr=nrow(gold_standard_network)
  nc=ncol(gold_standard_network)
  list_genes=union(gold_standard_network$V1,gold_standard_network$V2)
  numbergenes=length(list_genes)
  if (potentialtf=="known")
  {
	print("here we selected to work with known genes")
    mattf=as.matrix(gold_standard_network$V1[which(gold_standard_network$V3==1)])
    listtf<-unique(mattf)
    namefileglobalres=paste0("../res_in_sillico/insilico_size",sizenetwork,"_",network,"/globares_new_nbbootstrap_",nbBoobstrap,"_alphaenet_",tostringalphaenet,"_KOzcorenew_exponent_",tostringexponent,"_",sizenetwork,"_",network,".tsv")
  }
  else
  {
    listtf<- list_genes
    namefileglobalres=paste0("../res_in_sillico/insilico_size",sizenetwork,"_",network,"/globares_new_nbbootstrap_allgeneastf_",nbBoobstrap,"_alphaenet_",tostringalphaenet,"_KOzcorenew_exponent_",tostringexponent,"_",sizenetwork,"_",network,".tsv")
  }
  #tf=as.matrix(unique(gold_standard_network[gold_standard_network$V3==1,]$V1))
  print(paste("the number of transcription factor", length(listtf),sep = " "))

  matweight=matrix(0,nr = numbergenes, nc = numbergenes)
  rownames(matweight)<-list_genes
  colnames(matweight)<-list_genes
  for (g1 in list_genes)
  {
    for (g2 in list_genes)
    {
      #pval=matpval[g1,g2]
      #     # print(pval)
      # prob=ifelse(g1==g2,0,(intergral(beta=beta,pval=pval)))
      # matweight[g1,g2]=ifelse(prob>=0.8,0.0001,1)
      #matweight[g1,g2]=ifelse(g1==g2,0,1/#(integral(beta=beta,pval=pval,lambdamin=lambdamin,lambdamax=lambdamax)))
      matweight[g1,g2]=ifelse(g1==g2,0,1/abs(zscore[g1,g2])**exponent)
    }
  }
print(matweight)
  execTimeBenin<-system.time(edgeList<-applybootstrapbenin(X,nbBoobstrap=nbBoobstrap,matweightpk=matweight,sizenetwork = sizenetwork,
                                                           listtf = listtf,nbfolds = nbfolds,alphaenet = alphaenet,lmean=lmean,allgenes = list_genes))[[3]]
  print(paste("Benin took:", execTimeBenin," to proceed a network of size ", sizenetwork,sep = " "))
  tostringalphaenet=toString(alphaenet)
  tostringalphaenet=gsub(".","",tostringalphaenet,fixed=T)
  #namefileglobalres=paste0("../res_in_sillico/insilico_size",sizenetwork,"_",network,"/globares_new_nbbootstrap_",nbBoobstrap,"_alphaenet_",tostringalphaenet,"_KOzcore_exponent_",tostringexponent,"_",sizenetwork,"_",network,".tsv")
  global_res<-list2df(edgeList)
  global_res$TF<-as.character(global_res$TF)
  global_res$TG<-as.character(global_res$TG)
  test=global_res
  global_res[global_res[,1]==global_res[,2],3]=0
  global_res=subset(global_res,TF!=TG)
  global_res=global_res[order(global_res$W,decreasing = T),]
  savedata(data=global_res,file=namefileglobalres,sep="\t", colnames=FALSE, rownames=FALSE,quote=FALSE)
}
  
