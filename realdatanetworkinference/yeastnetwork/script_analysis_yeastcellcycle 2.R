library(org.Sc.sgd.db)
library(tidyr)
library(mice)
library(DMwR)
library(reshape2)
library(precrec)
library('foreach')
library('snow')
library('iterators')
library('doSNOW')
library('glmnet')
library('parallel')
library('boot')
library(doRNG)
library(dplyr)
### Function definition
savedata<- function(data,file, colnames=TRUE, rownames=TRUE,sep="\t",quote=FALSE,append=FALSE)
{
  write.table(data,file=file,sep=sep,append = append
              ,col.names=colnames,row.names = rownames,quote=FALSE)
}
proba<- function(beta,pval){
  proba2<- function(lambda){
    (lambda*exp(-lambda*pval)*beta)/((lambda*exp(-lambda*pval)*beta)+
                                       ((1-exp(-lambda))*(1-beta)))
  }
  return(proba2)
}
integral<- function(beta,pval,lambdamin,lambdamax){
  (integrate(proba(beta,pval),lower =lambdamin,upper = lambdamax)$value)/(lambdamax-lambdamin)
}
######### Beginning of the script

Sgdgenename <- org.Sc.sgdGENENAME
# Get the gene names that are mapped to an ORF identifier
mapped_genes <- mappedkeys(Sgdgenename)
# Convert to a list
list_mapped_genes <- as.list(Sgdgenename [mapped_genes])
### We match the set of gene that have previously been identified as cell cycle genes
rawtgs<-c("ALG7","CDC20","CDC21","CDC5","CDC6","CLB2","CLB5","CLN1","CLN2","CTS1","EGT2","FAR1","HTA1","PCL2","SIC1")
rawtfs<-c("ACE2","ASH1","FKH1","MBP1","MCM1","NDD1","STB1","SWI4","SWI5","SWI6")
rawallgenes<-c(rawtgs,rawtfs)
match_tg<-list_mapped_genes[which(list_mapped_genes%in%c("ALG7","CDC20","CDC21","CDC5","CDC6","CLB2","CLB5","CLN1","CLN2","CTS1","EGT2","FAR1","HTA1","PCL2","SIC1"))]
listtg<-names(match_tg)

##### We match the set of transcription factor that have previously been identified as cell cycle genes
match_tf<-list_mapped_genes[which(list_mapped_genes%in%c("ACE2","ASH1","FKH1","MBP1","MCM1","NDD1","STB1","SWI4","SWI5","SWI6"))]
#listtf<-names(match_tf)

allgenes_matching<-c(match_tg,match_tf)
##### we download the expression data matrix
listtfnames<-c("ACE2","ASH1","FKH1","MBP1","MCM1","NDD1","STB1","SWI4","SWI5","SWI6")
listtf<-unlist(lapply(listtfnames,FUN = function(x,db)names(db[(db==x)]),list_mapped_genes))
allgenes<-c(listtg,listtf)
ysexprdata<- read.csv("../data/time_expression_data/combined.txt",sep = "\t")
subysexprdata<-ysexprdata[(ysexprdata$X %in% allgenes),]
filteredexprdata<-subysexprdata[sapply(subysexprdata, function(x) !all(is.na(x)))] 
row.names(filteredexprdata )<-filteredexprdata [,1]
filteredexprdata<-filteredexprdata [,-1]
filteredexprdata<-t(filteredexprdata)
# imputed_filteredexprdata <- mice(filteredexprdata, m=5, maxit = 50, method = 'pmm', seed = 500,printFlag = FALSE)
# summary(imputed_filteredexprdata)
# finalexprdata<-complete(imputed_filteredexprdata,2)
cleanexprdata<-knnImputation(filteredexprdata,k=5)
cleanexprdata<-scale(x=cleanexprdata,scale = apply(cleanexprdata, 2, sd, na.rm = TRUE))
####### Reading location data 
locationdata<-read.csv("../data/location_data_saccarhomyces/binding_by_gene.tsv",sep = "\t")
filtered_locationdata<-locationdata[,-which(isNA(locationdata[1,]))]
filtered_locationdata<-filtered_locationdata[,-c(2:4)]
row.names(filtered_locationdata)<-filtered_locationdata[,1]
colnames(filtered_locationdata)<-as.matrix(filtered_locationdata[1,])
filtered_locationdata<-as.matrix(filtered_locationdata[,-1])
filtered_locationdata<-as.matrix(filtered_locationdata[-1,])
filtered_locationdata<-as.matrix(filtered_locationdata[allgenes,])
filtered_locationdata<-as.matrix(filtered_locationdata[,listtfnames])
colnames(filtered_locationdata)<-listtf
filtered_locationdata<-as.matrix(filtered_locationdata)
class(filtered_locationdata)<-"numeric"

# Transforming location data into probabilities
lambdamin=1
lambdamax=1000
beta=0.5
numbergenes<-length(allgenes)
nbtf<-length(listtf)
matweight=matrix(0,nr = numbergenes, nc = numbergenes)
rownames(matweight)<-allgenes
colnames(matweight)<-allgenes

for (g1 in allgenes)
{
  for (g2 in listtf)
  {
    pval=filtered_locationdata[g1,g2]
    #print(pval)
    # prob=ifelse(g1==g2,0,(intergral(beta=beta,pval=pval)))
    # matweight[g1,g2]=ifelse(prob>=0.8,0.0001,1)
    matweight[g1,g2]=ifelse(g1==g2,0,1/(integral(beta=beta,pval=pval,lambdamin=lambdamin,lambdamax=lambdamax)))
  }
}
######### Reading the gold standard network
gold_standard_network<-read.csv2("../data/gold_standard_ystract/RegulationTwoColumnTable_Documented_2013927.tsv",header = FALSE)
nbreglink<-dim(gold_standard_network)[1]
gold_standard_network<-cbind(gold_standard_network,rep(1,nbreglink))
colnames(gold_standard_network)<-c("TF","TG","W")
#### Transform into adjacency matrix
#adj_mat_gs_network<-gold_standard_network%>%dcast(TF~TG)
#row.names(adj_mat_gs_network)<-adj_mat_gs_network[,1]
#adj_mat_gs_network<-adj_mat_gs_network[,-1]
##### Setting parameters for BENIN

######## IMPORTANT IF IT DOES NOT WORK CHECK
#unlist(lapply(as.character(global_res[,1]),FUN=function(x,mapping)unlist(mapping[x],use.names = F),allgenes_matching),use.names = F)
nbBoobstrap=1000
alphaenet=0.4
sizenetwork=0
lmean=10
nbfolds=10
#normalize=TRUE,nbfolds,alphaenet,lmean=10,seed=500,listtf="",allgenes=""
tostringalphaenet<-toString(alphaenet)
tostringalphaenet=gsub(".","",tostringalphaenet,fixed=T)
edgeList<-bootstrapbenin(cleanexprdata,nbBoobstrap=nbBoobstrap,matweightpk=matweight,sizenetwork = sizenetwork,listtf = listtf,nbfolds = nbfolds,alphaenet = alphaenet,lmean=lmean,allgenes = allgenes)

global_res<-list2df(edgeList)
test<-global_res
#colnames(global_res)=c("TF","TG","W")
global_res[global_res$TF==global_res$TG,"W"]=0
#global_res[global_res[,1]==global_res[,2],3]=0
global_res<-global_res[order(global_res$W,decreasing = TRUE),]
#global_res<-global_res[which(global_res$PostProb>0),]
global_res$TF<-as.character(global_res$TF)
global_res$TG<-as.character(global_res$TG)
global_res[(global_res[,3]>=0.5),3]=1
global_res[(global_res[,3]<0.5),3]=0
#selctednetwork<-subset(global_res,W>=0.5)
# namefileglobalres<-paste("../res_in_sillico/","global_res_benin_yeast_alphaenet",tostringalphaenet
#                          ,'_nbbootstrap_',nbBoobstrap,'_lmean_',lmean,".tsv",sep = "")
# savedata(data=global_res,file=namefileglobalres,sep="\t",colnames=NA)
# res_for_eval=data.frame()
# res_for_eval=global_res
# res_for_eval[(res_for_eval[,3]>=0.5),3]=1
# res_for_eval[(res_for_eval[,3]<0.5),3]=0
# sgdALIAS <- org.Sc.sgdALIAS
# # Get the probe identifiers that are mapped to alias names
# mapped_probes <- mappedkeys(sgdALIAS)
# list_probes_names<-as.list(sgdALIAS[mapped_probes])
# selected_res_for_eval<-as.data.frame(res_for_eval[(res_for_eval[,3]==1),])
mappedinfreg<-unlist(allgenes_matching[as.character(global_res[,1])], use.names = FALSE)
mappedinftargen<-unlist(allgenes_matching[as.character(global_res[,2])], use.names = FALSE)
selected_res_for_eval_w_name<-global_res
#[(global_res2[,3]>=0.5),]
selected_res_for_eval_w_name$TF<-mappedinfreg
selected_res_for_eval_w_name$TG<-mappedinftargen
#colnames(selected_res_for_eval_w_name)<-colnames(gold_standard_network)
selected_res_for_eval_w_name$TF<-as.character(selected_res_for_eval_w_name$TF)
selected_res_for_eval_w_name$TG<-as.character(selected_res_for_eval_w_name$TG)
selected_res_for_eval_w_name<-subset(selected_res_for_eval_w_name,TG!=TF)
#### Subset of the gold standard network in cell cycle
sub_gold_standard_network=subset(gold_standard_network,TF%in%rawtfs& TG%in%rawallgenes)
sub_gold_standard_network$TF<-as.character(sub_gold_standard_network$TF)
sub_gold_standard_network$TG<-as.character(sub_gold_standard_network$TG)
sub_gold_standard_network<-subset(sub_gold_standard_network,TF!=TG)
#result<-left_join(selected_res_for_eval_w_name, gold_standard_network, by=c("TF","TG"),type="left",match="all")
result<-full_join(selected_res_for_eval_w_name,sub_gold_standard_network, by=c("TF","TG"))
result[is.na(result)]<-0
# pred<-prediction(result$W.x, result$W.y, label.ordering = NULL)
# curve <- performance(pred, "prec", "rec")
# namefilerecNetwork<-paste("../res_in_sillico/","global_res_for_eval",sizenetwork,"_",network,'_',file,"_benin_yeast_alphaenet",tostringalphaenet
#                           ,'_nbbootstrap_',nbBoobstrap,'_lmean_',lmean,'_',sizenetwork,"_",
#                           network,".tsv",sep = "")
# savedata(data=selected_res_for_eval,file=namefilerecNetwork, colnames=FALSE, rownames=FALSE,sep="\t",quote=FALSE)
# output<-res_for_eval%>%dcast(Regulator~TargetGene)
# output[is.na(output)]<-0
# row.names(output)<-output[,1]
# output<-output[,-1]


curve<-precrec::evalmod(scores=result$W.x, labels=result$W.y, mode="prcroc")
aucs <- precrec::auc(curve)
aucs$aucs
# source("https://bioconductor.org/biocLite.R")
# biocLite("limma")

# Location data alone
res_regnet_from_locdata=melt(filtered_locationdata, varname= c("TG","TF"),value.name = "W")
res_regnet_from_locdata=res_regnet_from_locdata[,c("TF","TG","W")]
res_regnet_from_locdata<-res_regnet_from_locdata[order(res_regnet_from_locdata$W,decreasing = FALSE),]
#res_regnet_from_locdata<-res_regnet_from_locdata[res_regnet_from_locdata$W<0.05,]
#res_regnet_from_locdata<-subset(res_regnet_from_locdata,TG%in%allgenes)
res_regnet_from_locdata[res_regnet_from_locdata$W>0.05,3]<-2
res_regnet_from_locdata[res_regnet_from_locdata$W<=0.05,3]<-1
res_regnet_from_locdata[res_regnet_from_locdata$W==2,3]<-0
#res_regnet_from_locdata<-res_regnet_from_locdata[res_regnet_from_locdata$W<=0.05,]
res_regnet_from_locdata$TF<-as.character(res_regnet_from_locdata$TF)
res_regnet_from_locdata$TG<-as.character(res_regnet_from_locdata$TG)
res_regnet_from_locdata<-subset(res_regnet_from_locdata,TF!=TG)

#res_regnet_from_locdata$W=1
mappedinfreg<-unlist(allgenes_matching[res_regnet_from_locdata[,1]], use.names = FALSE)
mappedinftargen<-unlist(allgenes_matching[res_regnet_from_locdata[,2]], use.names = FALSE)
colnames(res_regnet_from_locdata)<-colnames(gold_standard_network)
res_regnet_from_locdata$TF<-mappedinfreg
res_regnet_from_locdata$TG<-mappedinftargen
sub_gold_standard_network<-subset(gold_standard_network,TF%in%rawtfs& TG%in%rawallgenes)
sub_gold_standard_network$TF<-as.character(sub_gold_standard_network$TF)
sub_gold_standard_network$TG<-as.character(sub_gold_standard_network$TG)
sub_gold_standard_network<-subset(sub_gold_standard_network,TF!=TG)
result_loc_data<-full_join(res_regnet_from_locdata,sub_gold_standard_network, by=c("TF","TG"))
result_loc_data[is.na(result_loc_data)]<-0
curve_loc_data<-precrec::evalmod(scores=result_loc_data$W.x, labels=result_loc_data$W.y, mode="prcroc")
aucs_loc_data <- precrec::auc(curve_loc_data)
aucs_loc_data$aucs

global_resno_prior<-list2df(edgeListno_prior)
testno_prior<-global_resno_prior
colnames(global_resno_prior)=c("TF","TG","W")
global_resno_prior[(global_resno_prior[,3]>=0.5),3]<-1
global_resno_prior[(global_resno_prior[,3]<0.5),3]<-0
global_resno_prior[global_resno_prior[,1]==global_resno_prior[,2],3]<-0
mat_global_resnoprior<-global_resno_prior%>%dcast(TG~TF)
rownames(mat_global_resnoprior)<-mat_global_resnoprior[,1]
mat_global_resnoprior<-mat_global_resnoprior[,-1]
mat_global_resnoprior[ is.na(mat_global_resnoprior)]<-0
resfuncanalysis_noprior<-functionalAnalisisRandNetwork(allgenes,global_resno_prior)

global_res<-list2df(edgeList)
test<-global_res
colnames(global_res)=c("TF","TG","W")
global_res[(global_res[,3]>=0.5),3]<-1
global_res[(global_res[,3]<0.5),3]<-0
global_res[global_res[,1]==global_res[,2],3]<-0
mat_global_res<-global_res%>%dcast(TG~TF)
rownames(mat_global_res)<-mat_global_res[,1]
mat_global_res<-mat_global_res[,-1]
mat_global_res[ is.na(mat_global_res)]<-0
resfuncanalysis<-functionalAnalisisRandNetwork(allgenes,mat_global_res)




####### No prior 
edgeListno_prior<-bootstrapbenin(cleanexprdata,nbBoobstrap=nbBoobstrap,sizenetwork = sizenetwork,listtf = listtf,nbfolds = nbfolds,alphaenet = alphaenet,lmean=lmean,allgenes = allgenes)

global_resno_prior<-list2df(edgeListno_prior)
test<-global_res
#colnames(global_res)=c("TF","TG","W")
global_resno_prior[global_resno_prior$TF==global_resno_prior$TG,"W"]=0
global_resno_prior[global_resno_prior[,1]==global_resno_prior[,2],3]=0
global_resno_prior=global_resno_prior[order(global_resno_prior$W,decreasing = TRUE),]
#global_res<-global_res[which(global_res$PostProb>0),]

global_resno_prior[(global_resno_prior[,3]>=0.5),3]=1
global_resno_prior[(global_resno_prior[,3]<0.5),3]=0
#selctednetworkno_prior<-subset(global_resno_prior,W>=0.5)
# namefileglobalres<-paste("../res_in_sillico/","global_res_benin_yeast_alphaenet",tostringalphaenet
#                          ,'_nbbootstrap_',nbBoobstrap,'_lmean_',lmean,".tsv",sep = "")
# savedata(data=global_res,file=namefileglobalres,sep="\t",colnames=NA)
# res_for_eval=data.frame()
# res_for_eval=global_res
# res_for_eval[(res_for_eval[,3]>=0.5),3]=1
# res_for_eval[(res_for_eval[,3]<0.5),3]=0
# sgdALIAS <- org.Sc.sgdALIAS
# # Get the probe identifiers that are mapped to alias names
# mapped_probes <- mappedkeys(sgdALIAS)
# list_probes_names<-as.list(sgdALIAS[mapped_probes])
# selected_res_for_eval<-as.data.frame(res_for_eval[(res_for_eval[,3]==1),])
mappedinfregno_prior<-unlist(allgenes_matching[as.character(global_resno_prior[,1])], use.names = FALSE)
mappedinftargenno_prior<-unlist(allgenes_matching[as.character(global_resno_prior[,2])], use.names = FALSE)
selected_res_for_eval_w_nameno_prior<-global_resno_prior
#[(global_res2[,3]>=0.5),]
selected_res_for_eval_w_nameno_prior$TF<-mappedinfregno_prior
selected_res_for_eval_w_nameno_prior$TG<-mappedinftargenno_prior
#colnames(selected_res_for_eval_w_name)<-colnames(gold_standard_network)
selected_res_for_eval_w_nameno_prior$TF<-as.character(selected_res_for_eval_w_nameno_prior$TF)
selected_res_for_eval_w_nameno_prior$TG<-as.character(selected_res_for_eval_w_nameno_prior$TG)
selected_res_for_eval_w_nameno_prior<-subset(selected_res_for_eval_w_nameno_prior,TG!=TF)
#### Subset of the gold standard network in cell cycle
sub_gold_standard_network=subset(gold_standard_network,TF%in%rawtfs& TG%in%rawallgenes)
sub_gold_standard_network$TF<-as.character(sub_gold_standard_network$TF)
sub_gold_standard_network$TG<-as.character(sub_gold_standard_network$TG)
sub_gold_standard_network<-subset(sub_gold_standard_network,TF!=TG)
#result<-left_join(selected_res_for_eval_w_name, gold_standard_network, by=c("TF","TG"),type="left",match="all")
resultno_prior<-full_join(selected_res_for_eval_w_nameno_prior,sub_gold_standard_network, by=c("TF","TG"))
resultno_prior[is.na(resultno_prior)]<-0
# pred<-prediction(result$W.x, result$W.y, label.ordering = NULL)
# curve <- performance(pred, "prec", "rec")
# namefilerecNetwork<-paste("../res_in_sillico/","global_res_for_eval",sizenetwork,"_",network,'_',file,"_benin_yeast_alphaenet",tostringalphaenet
#                           ,'_nbbootstrap_',nbBoobstrap,'_lmean_',lmean,'_',sizenetwork,"_",
#                           network,".tsv",sep = "")
# savedata(data=selected_res_for_eval,file=namefilerecNetwork, colnames=FALSE, rownames=FALSE,sep="\t",quote=FALSE)
# output<-res_for_eval%>%dcast(Regulator~TargetGene)
# output[is.na(output)]<-0
# row.names(output)<-output[,1]
# output<-output[,-1]


curveno_prior<-precrec::evalmod(scores=resultno_prior$W.x, labels=resultno_prior$W.y, mode="prcroc")
aucsno_prior <- precrec::auc(curveno_prior)
aucsno_prior$aucs