library(PRROC)
library(precrec)
library(data.table)
library(dplyr)
library(tidyverse)
library(MRIaggr)
library(biomaRt)
library(org.Mm.eg.db)
library(igraph)
library(intergraph)
library(ggnet)
library(GGally)
library(ROCR)
library(biomaRt)
library(readr)
library(nVennR)
#
#home="/Volumes/DISK_IMG/GRN_inference_temp/src"
#home="/Volumes/Seagate\ Backup\ Plus\ Drive/PhD/GRN_inference_prime/src"
#home="/Users/Stephanie/Downloads/GRN_inference_temp"
home= "/Volumes/Seagate\ Backup\ Plus\ Drive/PhD/BENIN_git/"
setwd(home)
source("R/MRIaggr/Functions_Valid.R")
source("R/utile.R")
source("R/benin.R")
source("src/bootpracs.q")
source("src/bootfuns.q")

#loding goldstandard network
namefilegs="../data/data_human/final_data_hum_reg_network/Hela_data/final_human_goldstandard.txt"
#gold_standard_network<-fread(namefilegs,sep="\t", header=T, stringsAsFactors = F, quote = "",
#                             colClasses = c("character","character","numeric"))
#gold_standard_network<-
load(file="../data/data_human/final_data_hum_reg_network/Hela_data/final_human_goldstandard_improved_corrected.RData")
gold_standard_network<-final_human_gs
file_cell_cycle_genes<-"../data/data_human/final_data_hum_reg_network/Hela_data/considered_genes_hela_cell_cycle.txt"
cell_cycles_genes<-(read.table(file_cell_cycle_genes,header = F,quote="",sep="\t",stringsAsFactors = F))$V1
file_potential_tf<-"../data/data_human/final_data_hum_reg_network/Hela_data/Hela_cell_cycle_TF.txt"
potentialTFname<-as.vector(read.table(file_potential_tf,stringsAsFactors = F,quote = "",header = F)$V1)


#merging res from different cluster
#'functional',"KD","TFBS","Chipseq"
listpriordata= c("functional","KD","TFBS","Chipseq")
#listpriordata= "Chipseq"
#listpriordata="none"
nbcluster=5
alphaenet=0.9
exponent=1.5
nbBootstrap=5000
combmethod="max"
tostringalphaenet=toString(alphaenet)
tostringalphaenet=gsub(".","",tostringalphaenet,fixed=T)
tostringexponent=toString(exponent)
tostringexponent=gsub(".","",tostringexponent,fixed=T)
#allgenesinclusters<-list()
#reshelanetwork_with_TFBS_nbbootstrap_2500_alphaenet_09_exponent_09_helanetwork_1.txt
allfinalinferrednetwork<-data.frame()
prev=listpriordata[1]
for (priordata in listpriordata)
{
  finalinferrednetwork=data.frame()
  
  for (clusterid in 1:nbcluster)
  {
    fileresfile=paste0("../res_in_sillico/res_human/res_hela_network/cluster_",clusterid,"/reshelanetwork_with_",priordata,
                       "_nbbootstrap_",nbBootstrap,"_alphaenet_",tostringalphaenet,"_exponent_",tostringexponent,"_helanetwork_"
                       ,clusterid,".txt")
    finalres<-fread(fileresfile,sep="\t", header=T, stringsAsFactors = F, quote = "",
                    colClasses = c("character","character","numeric"))
    finalinferrednetwork<-rbind(finalinferrednetwork,finalres)
    #filecluster=paste0("../data/data_human/final_data_hum_reg_network/Hela_data/modules/cluster_",clusterid,"_comp_clust.txt")
    #filecluster=paste0("../data/data_human/final_data_hum_reg_network/Hela_data/modules/",clusterid,".txt")
    #genesincluster<-as.vector(read.table(filecluster,stringsAsFactors = F,header = F,quote="")$V1)
    #allgenesinclusters<-c(allgenesinclusters,genesincluster)
    #print(setdiff(potentialTFname,finalinferrednetwork$TF))
  }
  if (dim(allfinalinferrednetwork)[1]==0)
  {
    allfinalinferrednetwork<-finalinferrednetwork
  }
  else
  {
    allfinalinferrednetwork<-merge(allfinalinferrednetwork,finalinferrednetwork,by=c("TF","TG"),all = T,suffixes = c(prev,priordata))
  }
  prev=priordata
}

if(length(priordata)==1 & length(listpriordata)==1)
{
  colnames(allfinalinferrednetwork)<-c("TF","TG","infW")
}else
{
  if (combmethod=="max")
  {
    allfinalinferrednetwork$infW<-apply(allfinalinferrednetwork[,3:6],1,FUN=function(x){max(x,na.rm=T)})
  }else
  {
    if  (combmethod=="mean")
    {
      allfinalinferrednetwork$infW<-apply(allfinalinferrednetwork[,3:6],1,FUN=function(x){mean(x,na.rm=T)})
    }
  }
}
#finalinferrednetwork<-finalinferrednetwork[order(finalinferrednetwork$W,decreasing = T),]
#finalresforeval<-merge(gold_standard_network,finalinferrednetwork,by=c("TF","TG"),sort=F,all.x=T)
#finalresforeval[is.na(finalresforeval$W.y),"W.y"]<-0
#finalresforeval<-subset(finalresforeval,!is.na(W.x))

allfinalinferrednetwork<-allfinalinferrednetwork[order(allfinalinferrednetwork$infW,decreasing = T),]




uniquegold_standard_network<-gold_standard_network%>%distinct(TF,TG,.keep_all =T)
#newgold_standard_network=reshape2::dcast(uniquegold_standard_network,TF~TG,fill=0)
#uniquegold_standard_network=reshape2::melt(newgold_standard_network)
#colnames(uniquegold_standard_network)<-c("TF","TG","W")
#uniquegold_standard_network<-gold_standard_network%>%distinct(TF,TG,.keep_all =T)
#uniquegold_standard_network<-uniquegold_standard_network%>%distinct(TF,TG,.keep_all =T)
#uniquegold_standard_network<-subset(uniquegold_standard_network,TF!=TG)
finalresforeval<-merge(uniquegold_standard_network,allfinalinferrednetwork[,c("TF","TG","infW")],by=c("TF","TG"),sort=F,all.x=T)
finalresforeval<-subset(finalresforeval,TF!=TG)
finalresforeval[is.na(finalresforeval$infW),"infW"]<-0

#finalresforeval<-subset(finalresforeval,!is.na(maxcomb))
#f
otherfinalresforeval<-merge(uniquegold_standard_network,allfinalinferrednetwork[,c("TF","TG","infW")],
                            by=c("TF","TG"),sort=F,all=T)

subotherfinalresforeval<-subset(otherfinalresforeval,is.na(W)&TF%in%potentialTFname)
thres<-0.5
finalresforeval[finalresforeval$infW<thres,"infW"]<-0


allfinalinferrednetwork_foranalysis<-subset(allfinalinferrednetwork,infW>=thres)
########## To register
fileresbenin=paste0("../res_in_sillico/res_human/res_hela_network/Hela_network_all_prior_BENIN_alphaenet",
                    tostringalphaenet,"_exponent_",tostringexponent,"_nbbootstrap_",nbBootstrap,"_allinfo.txt")

reg_otherfinalresforeval<-merge(uniquegold_standard_network,allfinalinferrednetwork,
                            by=c("TF","TG"),sort=F,all=T)
reg_otherfinalresforeval<-subset(reg_otherfinalresforeval,TF%in%potentialTFname)

#write.table(reg_otherfinalresforeval,fileresbenin,col.names = T,row.names = F, quote = F,sep = "\t")
# finalresforeval[,"W"]<-1
# finalresforeval[2:10,"W"]<-0
# finalresforeval[,"infW"]<-0
AUPRC=MLmetrics::PRAUC(finalresforeval$infW,finalresforeval$W)
print(paste0("The AUPR score for expression alone is ",AUPRC," we use combination with ",combmethod))
#print(AUPRC)
AUROC=MLmetrics::AUC(finalresforeval$infW,finalresforeval$W)
print(paste0("The AUROC score for expression alone is ",AUROC," we use combination with ",combmethod))
#print(AUROC)
pred <- prediction(finalresforeval$infW,finalresforeval$W)
roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")


####### For roc curve
if (length(listpriordata)==1)
{
  title=paste0("ROCR TS + ",listpriordata[1])
}else{
  title=paste0("ROCR TS + KD + Chipseq+ Functional + TFBS")
}
plot(roc.perf,main=title,colorize=T)
abline(a=0, b= 1)
auc.perf = performance(pred, measure = "auc")
auc.perf@y.values
pr.perf = performance(pred, measure="prec", x.measure="rec")

##### For PR curve
if (length(listpriordata)==1)
{
  title=paste0("PRCR TS + ",listpriordata[1])
}else{
  title=paste0("PRCR TS + KD + Chipseq+ Functional + TFBS")
}
plot(pr.perf,main=title,colorize=T)
score=modcalcAUPRC(performance = pr.perf, subdivisions = 100000,method = "integrate")
print(score)
x<-pr.perf@x.values[[1]]
x<-x[2:length(x)]
y<-pr.perf@y.values[[1]]
y<-y[2:length(y)]
pracma::trapz(x,y)

pROC_obj <-pROC::roc(finalresforeval$W,finalresforeval$infW,smoothed = TRUE,
                     # arguments for ci
                     ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                     # arguments for plot
                     plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                     print.auc=TRUE, show.thres=TRUE)


fg <- subset(finalresforeval,W==1)$infW
bg <- subset(finalresforeval,W==0)$infW
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
plot(pr)

roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
plot(roc)

#perfXY <- ROCR::performance(ROCR::prediction(finalresforeval$infW,finalresforeval$W), x.measure = "rec", measure = "prec")
#modcalcAUPRC(performance = perfXY, subdivisions = 10000)




################################# Orthology network


################ Getting all TF from human
humantffile<-"../data/data_human/final_data_hum_reg_network/Hela_data/human_TF_list.txt"
humantfinfo<-read.table(humantffile,stringsAsFactors = F,header = T, quote="",sep="\t",fill=T)
sub_humantfinfo<-subset(humantfinfo,Is_TF=="Yes")
humantffile2<-"../data/data_human/final_data_hum_reg_network/Hela_data/Homo_sapiens_TF_animalTFDB.txt"
humantfinfo2<-read.table(humantffile2,stringsAsFactors = F,header = T, quote="",sep="\t",fill=T)
allHumanTF<-union(sub_humantfinfo$Name,humantfinfo2$Symbol)
 ##### Getting the orthologs genes
selectrowmultannot<-function(x,sep=",")
{
  listannot<-unlist(strsplit(x[4],split=sep,fixed = T)[[1]])
  nbannot<-length(listannot)
#print(x)
  if (nbannot>1)
  {
    #return (list(x[1],x[2],x[3],x[4]))
    return (as.vector(c(x[1],x[2],x[3],x[4],x[5])))
  }
}

splitrow<-function(x,data )
{
  nbrow<-nrow(data)
  pos=nbrow+1
  listannot<-unlist(strsplit(x[3],",",fixed = T)[[1]])
  nbannot<-length(listannot)

  if (nbannot>1)
  {
    query<-x[1]
    species<-x[2]
    genename<-x[4]
    for (indx in 1:nbannot)
    {
      ortholog<-listannot[indx]
      #print(ortholog)
      #data[pos,]<-c(query,species,ortholog,genename)
      #pos=pos+1
    }
  }
}
filemappinguprotgenename<-"../res_in_sillico/human_orth_eggnog/mapping_uniprot_genename.txt"
mappinguprotgenename<-read.table(filemappinguprotgenename,stringsAsFactors = F,sep="\t",quote = "",header = T)
colnames(mappinguprotgenename)<-c("uniprotID","genename")
fileortholog<-"../res_in_sillico/human_orth_eggnog/hum-cellcycle-ouput.emapper.predict_orthologs"
orthologdata<-read.table(fileortholog,header = T, stringsAsFactors = F, quote="\"'", sep="\t")
#genename<-unlist(lapply(orthologdata$Query,FUN=function(x){strsplit(x,"|",fixed = T)[[1]][3]}))
#genename<-unlist(lapply(genename,FUN=function(x){strsplit(x,"_",fixed = T)[[1]][1]}))
#orthologdata$genename<-genename
orthologdata$uniprotID<-unlist(lapply(orthologdata$Query,FUN=function(x){strsplit(x,"|",fixed = T)[[1]][2]}))
orthologdata<-merge(orthologdata,mappinguprotgenename,by="uniprotID",all.x=T)
# ############################################################Zebra
# ############ Orthology from zebrafish
# orthologdata_from_zebra<-subset(orthologdata,Species%in%"7955")
# #genename<-unlist(lapply(orthologdata_from_zebra$Query,FUN=function(x){strsplit(x,"|",fixed = T)[[1]][3]}))
# #genename<-unlist(lapply(genename,FUN=function(x){strsplit(x,"_",fixed = T)[[1]][1]}))
# #orthologdata_from_zebra$genename<-genename
# multannot<-apply(orthologdata_from_zebra,1,FUN=selectrowmultannot)
# multannot<-multannot[!unlist(lapply(multannot, is.null))]
# 
# for (annotinfo in multannot)
# {
#   nbrow<-nrow(orthologdata_from_zebra)
#   pos=nbrow+1
#   listannot<-unlist(strsplit(annotinfo[4],",",fixed = T)[[1]])
#   nbannot<-length(listannot)
#   
#   if (nbannot>1)
#   {
#     uniprotID<-annotinfo[1]
#     query<-annotinfo[2]
#     species<-annotinfo[3]
#     genename<-annotinfo[5]
#     for (indx in 1:nbannot)
#     {
#       ortholog<-listannot[indx]
#       #print(ortholog)
#       orthologdata_from_zebra[pos,]<-c( uniprotID,query,species,ortholog,genename)
#       #pos=pos+1
#     }
#   }
# }
# orthologdata_from_zebra<-orthologdata_from_zebra[!str_detect(orthologdata_from_zebra$Orthologs,","),]
# orthologdata_from_zebra$ID.ortholog<-unlist(lapply(orthologdata_from_zebra$Orthologs,FUN=function(x){strsplit(x,".",fixed = T)[[1]][2]}))
# 
# 
# 
# ############ Info TF zebra fish
# fileTFInfomotifzebra<-"../data/zebbrafish_info_regulation/Danio_rerio_TF.txt"
# TFInfomotifzebra<-read.table(fileTFInfomotifzebra,header = T,stringsAsFactors = F,sep = "\t",quote="")
# sub_TFInfomotifzebra<-TFInfomotifzebra[,c("Symbol","Ensembl","Family","Protein","Entrez.ID")]
# multannot<-apply(sub_TFInfomotifzebra,1,FUN=selectrowmultannot,sep=";")
# multannot<-multannot[!unlist(lapply(multannot, is.null))]
# genewithmultprot<-unlist(lapply(multannot,FUN=function(x)x[[1]]))
# clean_sub_TFInfomotifzebra<-subset(sub_TFInfomotifzebra,!Symbol%in%genewithmultprot)
# # 
# # 
# for (annotinfo in multannot)
# {
#   nbrow<-nrow(clean_sub_TFInfomotifzebra)
#   pos=nbrow+1
#   listannot<-unlist(strsplit(annotinfo[4],";",fixed = T)[[1]])
#   nbannot<-length(listannot)
#   if (nbannot>1)
#   {
#     
#     Symbol<-annotinfo[1]
#     Ensembl<-annotinfo[2]
#     Family<-annotinfo[3]
#     Entrez.ID<-annotinfo[5]
#     for (indx in 1:nbannot)
#     {
#       protein<-listannot[indx]
#       if (protein!="")
#       {
#         #print(c(Species,Symbol,protein,Entrez.ID))
#         clean_sub_TFInfomotifzebra[pos,]<-c(Symbol,Ensembl,Family,protein,Entrez.ID)
#       }
#       
#       
#       #print(ortholog)
#       
#       #pos=pos+1
#     }
#   }
# }
# clean_sub_TFInfomotifzebra$Protein<-gsub(";","",clean_sub_TFInfomotifzebra$Protein)
# sub_clean_sub_TFInfomotifzebra<-subset(clean_sub_TFInfomotifzebra,Protein%in%orthologdata_from_zebra$ID.ortholog)
# 
# ####### Reg network from zebra
# 
# filezebranetwork<-"../data/zebbrafish_info_regulation/7955.protein.actions.v11.0.txt"
# zebranetwork<-fread(filezebranetwork,header = T, stringsAsFactors = F,sep="\t",
#                         colClasses = c("character","character","character","character","character","character","numeric"))
# sub_zebranetwork<-zebranetwork%>%subset(mode%in%c("binding","expression") & score>=200 &item_id_a%in%orthologdata_from_zebra$Orthologs&
#                                               item_id_b%in%orthologdata_from_zebra$Orthologs&(is_directional%in%"t"|a_is_acting%in%"t"))
# 
# sub_zebranetwork<-sub_zebranetwork[,c("item_id_a","item_id_b")]
# colnames(sub_zebranetwork)<-c("TF","TG")
# 
# sub_zebranetwork$TF<-unlist(lapply(sub_zebranetwork$TF,FUN=function(x){strsplit(x,".",fixed = T)[[1]][2]}))
# sub_zebranetwork$TG<-unlist(lapply(sub_zebranetwork$TG,FUN=function(x){strsplit(x,".",fixed = T)[[1]][2]}))
# #sub_zebranetwork<-subset(sub_zebranetwork,TF%in%sub_clean_sub_TFInfomotifzebra$Protein|TG%in%sub_clean_sub_TFInfomotifzebra$Protein)
# ortho_human_network_from_zebra<-merge(sub_zebranetwork,orthologdata_from_zebra,
#                                       by.x="TF",by.y="ID.ortholog", all.x=T)
# 
# final_ortho_human_network_from_zebra<-merge(ortho_human_network_from_zebra,orthologdata_from_zebra,
#                                             by.x="TG",by.y="ID.ortholog", all.x=T,suffixes=c(".TF",".TG"))
# 
# final_ortho_human_network_from_zebra<-final_ortho_human_network_from_zebra%>%subset(!is.na(genename.TF)&!is.na(genename.TG))
# 
# final_ortho_human_network_from_zebra<-final_ortho_human_network_from_zebra[,c("genename.TF","genename.TG")]
# colnames(final_ortho_human_network_from_zebra)<-c("TF","TG")
# 
# #test<-sub_zebranetwork%>%subset(TF%in%sub_TFInfomotifzebra$Protein)
# sub_final_ortho_human_network_from_zebra<-subset(final_ortho_human_network_from_zebra,TF%in%intersect(humantfinfo$Name,cell_cycles_genes))
# sub_final_ortho_human_network_from_zebra<-sub_final_ortho_human_network_from_zebra%>%distinct(TF,TG,.keep_all =T)
# sub_final_ortho_human_network_from_zebra$infW<-1
# info_network_ortho_expr_zebra<-merge(allfinalinferrednetwork[,c("TF","TG","infW")],
#                                sub_final_ortho_human_network_from_zebra,by=c("TF","TG"),sort=F,all.y=T,suffixes=c(".expr",".orth"))
# info_network_ortho_expr_zebra<-merge(uniquegold_standard_network,info_network_ortho_expr_zebra,by=c("TF","TG"),sort=F,all.y=T)
# 
# 
# 
# #### definecolour for plotting
# color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
# ######choice from 334, 26,142,74,232,329,411
# 
# 
# ###### building network
# #network_gs<-graph_from_data_frame(info_network_ortho_expr[!is.na(W)&W>0,c("TF","TG","W")])
# network_gs<-graph_from_data_frame(subset(info_network_ortho_expr_zebra,!is.na(W)&W>0,select=c("TF","TG","W")),directed = T)
# #darkcyan
# network_gs<-network_gs%>%set_edge_attr("color", value = color[74])
# #graph_network_expr<-graph_from_data_frame(info_network_ortho_expr[!is.na(infW.expr)& infW.expr>=0.5,c("TF","TG","infW.expr")],directed = T)
# graph_network_expr<-graph_from_data_frame(subset(info_network_ortho_expr_zebra,!is.na(infW.expr)& infW.expr>=0.5
#                                                  ,select=c("TF","TG","infW.expr")),directed = T)
# 
# 
# #purple
# graph_network_expr<-graph_network_expr%>%set_edge_attr("color", value = color[329])
# #graph_network_ortho<-graph_from_data_frame(info_network_ortho_expr[!is.na(infW.orth),c("TF","TG","infW.orth")],directed = T)
# graph_network_ortho<-graph_from_data_frame(subset(info_network_ortho_expr_zebra,!is.na(infW.orth),select=c("TF","TG","infW.orth")),directed = T)
# 
# #green
# graph_network_ortho<-graph_network_ortho%>%set_edge_attr("color", value = color[142])
# allnetwork<-graph_network_ortho+graph_network_expr+network_gs
# nbedges<-length(E(allnetwork)$color)
# edgecol<-unlist(lapply(1:nbedges, FUN =getcolor,E(allnetwork)$color,E(allnetwork)$color_1,E(allnetwork)$color_2,color ))
# E(allnetwork)$color<-edgecol
# ### Saving the orthonetwork from mouse to a file
# model="zebra"
# namefilegraph=paste0("../res_in_sillico/res_orthology/res_from_orthology_with_",model,"_withexpr.gml")
# #write_graph(allnetwork,namefilegraph,"gml")
# 
# 
# 
# 
# 
# #####################################################################
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #####################################################################Saccharomyces
# 
# ############ Orthology from saccarhomyces
# orthologdata_from_saccarhomyce<-subset(orthologdata,Species%in%"4932")
# genename<-unlist(lapply(orthologdata_from_saccarhomyce$Query,FUN=function(x){strsplit(x,"|",fixed = T)[[1]][3]}))
# genename<-unlist(lapply(genename,FUN=function(x){strsplit(x,"_",fixed = T)[[1]][1]}))
# orthologdata_from_saccarhomyce$genename<-genename
# multannot<-apply(orthologdata_from_saccarhomyce,1,FUN=selectrowmultannot)
# multannot<-multannot[!unlist(lapply(multannot, is.null))]
# 
# 
# for (annotinfo in multannot)
# {
#   nbrow<-nrow(orthologdata_from_saccarhomyce)
#   pos=nbrow+1
#   listannot<-unlist(strsplit(annotinfo[3],",",fixed = T)[[1]])
#   nbannot<-length(listannot)
#   
#   if (nbannot>1)
#   {
#     query<-annotinfo[1]
#     species<-annotinfo[2]
#     genename<-annotinfo[4]
#     for (indx in 1:nbannot)
#     {
#       ortholog<-listannot[indx]
#       #print(ortholog)
#       orthologdata_from_saccarhomyce[pos,]<-c(query,species,ortholog,genename)
#       #pos=pos+1
#     }
#   }
# }
# orthologdata_from_saccarhomyce<-orthologdata_from_saccarhomyce[!str_detect(orthologdata_from_saccarhomyce$Orthologs,","),]
# orthologdata_from_saccarhomyce$ID.ortholog<-unlist(lapply(orthologdata_from_saccarhomyce$Orthologs,FUN=function(x){strsplit(x,".",fixed = T)[[1]][2]}))
# 
# ##### Regulatory network saccharomyces
# filesaccregnetwork<-"../data/sacc_info_regulation/RegulationTwoColumnTable_Documented_202066_1618_1097109570.tsv"
# saccregnetwork<-read.table(filesaccregnetwork,header = F, stringsAsFactors = F,sep=";")
# saccregnetwork$V1<-toupper(saccregnetwork$V1)
# saccregnetwork$V2<-toupper(saccregnetwork$V2)
# 
# 
# 
# #####################################################################
# 
# 



#####################################################################Mouse
############ Orthology from from mouse
orthologdata_from_mouse<-subset(orthologdata,Species%in%"10090")
#genename<-unlist(lapply(orthologdata_from_mouse$Query,FUN=function(x){strsplit(x,"|",fixed = T)[[1]][3]}))
#genename<-unlist(lapply(genename,FUN=function(x){strsplit(x,"_",fixed = T)[[1]][1]}))
#orthologdata_from_mouse$genename<-genename
multannot<-apply(orthologdata_from_mouse,1,FUN=selectrowmultannot)
multannot<-multannot[!unlist(lapply(multannot, is.null))]
#lapply(multannot,splitrow, data=orthologdata_from_mouse)

for (annotinfo in multannot)
{
  nbrow<-nrow(orthologdata_from_mouse)
  pos=nbrow+1
  listannot<-unlist(strsplit(annotinfo[3],",",fixed = T)[[1]])
  nbannot<-length(listannot)

  if (nbannot>1)
  {
    query<-annotinfo[1]
    species<-annotinfo[2]
    genename<-annotinfo[4]
    for (indx in 1:nbannot)
    {
      ortholog<-listannot[indx]
      #print(ortholog)
      orthologdata_from_mouse[pos,]<-c(query,species,ortholog,genename)
      #pos=pos+1
    }
  }
}
orthologdata_from_mouse<-orthologdata_from_mouse[!str_detect(orthologdata_from_mouse$Orthologs,","),]
orthologdata_from_mouse$ID.ortholog<-unlist(lapply(orthologdata_from_mouse$Orthologs,FUN=function(x){strsplit(x,".",fixed = T)[[1]][2]}))




# ############# List of tfs
filelisttfmouse<-"../data/mouse_info_regulation/Mus_musculus_TF.txt"
listtfmouseinfo<-read.table(filelisttfmouse,stringsAsFactors = F,header = T,sep="\t",quote = "")
listtfmouseinfo$Symbol<-toupper(listtfmouseinfo$Symbol)
listtfmouse<-listtfmouseinfo$Symbol

######### Reg network mouse

filemouseregnetwork<-"../data/mouse_info_regulation/trrust_rawdata.mouse.tsv"
mouseregnetwork<-read.table(filemouseregnetwork,header = F, stringsAsFactors = F,sep="\t")
mouseregnetwork$V1<-toupper(mouseregnetwork$V1)
mouseregnetwork$V2<-toupper(mouseregnetwork$V2)
sub_mouseregnetwork<-subset(mouseregnetwork,V1%in%orthologdata_from_mouse$genename &V2%in%orthologdata_from_mouse$genename)
sub_mouseregnetwork<-sub_mouseregnetwork[,c("V1","V2")]
colnames(sub_mouseregnetwork)<-c("TF","TG")
sub_mouseregnetwork<-subset(sub_mouseregnetwork,TF%in%listtfmouse)

filemouseregnetwork2<-"../data/mouse_info_regulation/mouse/mouse.source.txt"
mouseregnetwork2<-read.table(filemouseregnetwork2,header = F, stringsAsFactors = F,sep="\t")
mouseregnetwork2$V1<-toupper(mouseregnetwork2$V1)
mouseregnetwork2$V3<-toupper(mouseregnetwork2$V3)
sub_mouseregnetwork2<-subset(mouseregnetwork2,V1%in%orthologdata_from_mouse$genename &V3%in%orthologdata_from_mouse$genename)
sub_mouseregnetwork2<-sub_mouseregnetwork2[,c("V1","V3")]
colnames(sub_mouseregnetwork2)<-c("TF","TG")
sub_mouseregnetwork2<-subset(sub_mouseregnetwork2,TF%in%listtfmouse)



filemouseregnetwork3<-"../data/mouse_info_regulation/new_kegg.mouse.reg.direction_new.txt"
mouseregnetwork3<-read.table(filemouseregnetwork3,header = T, stringsAsFactors = F,sep="\t")
mouseregnetwork3$TF<-toupper(mouseregnetwork3$TF)
mouseregnetwork3$Target<-toupper(mouseregnetwork3$Target)
sub_mouseregnetwork3<-subset(mouseregnetwork3,TF%in%orthologdata_from_mouse$genename &Target%in%orthologdata_from_mouse$genename)
sub_mouseregnetwork3<-sub_mouseregnetwork3[,c("TF","Target")]
colnames(sub_mouseregnetwork3)<-c("TF","TG")
sub_mouseregnetwork3<-subset(sub_mouseregnetwork3,TF%in%listtfmouse)

filemouseregnetwork4<-"../data/mouse_info_regulation/10090.protein.actions.v11.0.txt"
mouseregnetwork4<-fread(filemouseregnetwork4,header = T, stringsAsFactors = F,sep="\t",
                        colClasses = c("character","character","character","character","character","character","numeric"))
mouseregnetwork4<-mouseregnetwork4%>%subset(mode%in%c("binding","expression") & score>=500 &item_id_a%in%orthologdata_from_mouse$Orthologs&
                                              item_id_b%in%orthologdata_from_mouse$Orthologs&is_directional%in%"t"&a_is_acting%in%"t")

sub_mouseregnetwork4<-mouseregnetwork4[,c("item_id_a","item_id_b")]
colnames(sub_mouseregnetwork4)<-c("TF","TG")
sub_mouseregnetwork4$TF<-unlist(lapply(sub_mouseregnetwork4$TF,FUN=function(x){strsplit(x,".",fixed = T)[[1]][2]}))
sub_mouseregnetwork4$TG<-unlist(lapply(sub_mouseregnetwork4$TG,FUN=function(x){strsplit(x,".",fixed = T)[[1]][2]}))
#sub_mouseregnetwork4<-subset(sub_mouseregnetwork4,TF%in%listtfmouse)

#mouseregnetwork3$TF<-toupper(mouseregnetwork3$TF)
#mouseregnetwork3$Target<-toupper(mouseregnetwork3$Target)
#sub_mouseregnetwork3<-subset(mouseregnetwork3,TF%in%orthologdata_from_mouse$genename &Target%in%orthologdata_from_mouse$genename)
#sub_mouseregnetwork3<-sub_mouseregnetwork3[,c("TF","Target")]
#colnames(sub_mouseregnetwork3)<-c("Gene1","Gene2")
final_mouseregulatorynetwork<-sub_mouseregnetwork
final_mouseregulatorynetwork<-rbind(final_mouseregulatorynetwork,sub_mouseregnetwork2)
final_mouseregulatorynetwork<-rbind(final_mouseregulatorynetwork,sub_mouseregnetwork3)
final_mouseregulatorynetwork[final_mouseregulatorynetwork$TG%in%"INADL","TG"]<-"PATJ"


mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL")
m_ensembl = useDataset(dataset = "mmusculus_gene_ensembl", mart = mart)
m_genes_mapping<-getBM(attributes = c("ensembl_peptide_id","wikigene_name"),filter="wikigene_name",
                       values=union(final_mouseregulatorynetwork$TF,final_mouseregulatorynetwork$TG),mart=m_ensembl)
m_genes_mapping$wikigene_name<-toupper(m_genes_mapping$wikigene_name)
m_genes_mapping<-subset(m_genes_mapping,ensembl_peptide_id!="")
final_mouseregulatorynetwork<-merge(final_mouseregulatorynetwork,m_genes_mapping,by.x="TF",
                                    by.y="wikigene_name",all.x=T)
final_mouseregulatorynetwork<-merge(final_mouseregulatorynetwork,m_genes_mapping,by.x="TG",
                                    by.y="wikigene_name",all.x=T,suffixes=c("TF","TG"))
final_mouseregulatorynetwork<-final_mouseregulatorynetwork[,c("ensembl_peptide_idTF","ensembl_peptide_idTG")]
colnames(final_mouseregulatorynetwork)<-c("TF","TG")
#h_ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
final_mouseregulatorynetwork<-rbind(final_mouseregulatorynetwork,sub_mouseregnetwork4)
final_mouseregulatorynetwork<-final_mouseregulatorynetwork%>%distinct(TF,TG,.keep_all = T)

mappinglisttfmouse<-getBM(attributes = c("ensembl_peptide_id","wikigene_name"),filter="wikigene_name",
                          values=listtfmouse,mart=m_ensembl)

subfinal_mouseregulatorynetwork<-subset(final_mouseregulatorynetwork,
                                        TF%in%mappinglisttfmouse$ensembl_peptide_id)
#subfinal_mouseregulatorynetwork<-final_mouseregulatorynetwork
ortho_human_network_from_mouse<-merge(subfinal_mouseregulatorynetwork,orthologdata_from_mouse,
                                      by.x="TF",by.y="ID.ortholog",all.x=T)
ortho_human_network_from_mouse<-subset(ortho_human_network_from_mouse,!is.na(genename))
#ortho_human_network_from_mouse<-merge(ortho_human_network_from_mouse,orthologdata_from_mouse,
#                                      by.x="TG",by.y="ID.ortholog", all.x=T, suffixes=c("TF","TG"))
final_ortho_human_network_from_mouse<-merge(ortho_human_network_from_mouse,orthologdata_from_mouse,
                                            by.x="TG",by.y="ID.ortholog", suffixes=c("TF","TG"),all.x=T)
final_ortho_human_network_from_mouse<-subset(final_ortho_human_network_from_mouse,!is.na(genenameTG)&!is.na(genenameTF))
final_ortho_human_network_from_mouse<-final_ortho_human_network_from_mouse[,c("genenameTF","genenameTG")]
colnames(final_ortho_human_network_from_mouse)<-c("TF","TG")
final_ortho_human_network_from_mouse<-final_ortho_human_network_from_mouse%>%distinct(TF,TG,.keep_all = T)
final_ortho_human_network_from_mouse$infW<-1
sub_final_ortho_human_network_from_mouse<-subset(final_ortho_human_network_from_mouse,TF%in%intersect(allHumanTF,cell_cycles_genes))
sub_final_ortho_human_network_from_mouse<-sub_final_ortho_human_network_from_mouse%>%distinct(TF,TG,.keep_all =T)


#finalresforeval_w_ortho<-rbind(allfinalinferrednetwork[,c("TF","TG","infW")],final_ortho_human_network_from_mouse)

#uniquegold_standard_network<-subset(uniquegold_standard_network,TF!=TG)
#finalresforeval_w_ortho<-merge(uniquegold_standard_network,finalresforeval_w_ortho,by=c("TF","TG"),sort=F,all.x=T)
info_network_ortho_expr<-merge(allfinalinferrednetwork_foranalysis[,c("TF","TG","infW")],
                               sub_final_ortho_human_network_from_mouse,by=c("TF","TG"),sort=F,all.y=T,suffixes=c(".expr",".orth"))
info_network_ortho_expr<-merge(uniquegold_standard_network,info_network_ortho_expr,by=c("TF","TG"),sort=F,all.y=T)


#### definecolour for plotting
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
######choice from 334, 26,142,74,232,329,411


###### building network
#network_gs<-graph_from_data_frame(info_network_ortho_expr[!is.na(W)&W>0,c("TF","TG","W")])
network_gs<-graph_from_data_frame(subset(info_network_ortho_expr,!is.na(W),select=c("TF","TG","W")),directed = T)
#darkcyan
network_gs<-network_gs%>%set_edge_attr("color", value = color[74])
#graph_network_expr<-graph_from_data_frame(info_network_ortho_expr[!is.na(infW.expr)& infW.expr>=0.5,c("TF","TG","infW.expr")],directed = T)
graph_network_expr<-graph_from_data_frame(subset(info_network_ortho_expr,!is.na(infW.expr)
                                                 ,select=c("TF","TG","infW.expr")),directed = T)


#purple
graph_network_expr<-graph_network_expr%>%set_edge_attr("color", value = color[329])
#graph_network_ortho<-graph_from_data_frame(info_network_ortho_expr[!is.na(infW.orth),c("TF","TG","infW.orth")],directed = T)
graph_network_ortho<-graph_from_data_frame(subset(info_network_ortho_expr,!is.na(infW.orth),select=c("TF","TG","infW.orth")),directed = T)

#green
graph_network_ortho<-graph_network_ortho%>%set_edge_attr("color", value = color[142])
allnetwork<-graph_network_ortho+graph_network_expr+network_gs
nbedges<-length(E(allnetwork)$color)
edgecol<-unlist(lapply(1:nbedges, FUN =getcolor,E(allnetwork)$color,E(allnetwork)$color_1,E(allnetwork)$color_2,color ))
E(allnetwork)$color<-edgecol
### Saving the orthonetwork from mouse to a file
model="mouse"
namefilegraph=paste0("../res_in_sillico/res_orthology/res_from_orthology_with_",model,"_finalaftercorrection_withexpr.gml")
write_graph(allnetwork,namefilegraph,"gml")





######### We want to evaluate
finalresforeval_ortho_expr<-merge(final_ortho_human_network_from_mouse,finalresforeval,
                                  by=c("TF","TG"),suffixes=c(".ortho",".expr"),all=T)

sub_finalresforeval_ortho_expr<-subset(finalresforeval_ortho_expr,!is.na(W))
combmethod="mean"
if (combmethod=="max")
{
  sub_finalresforeval_ortho_expr$infWcom<-apply(sub_finalresforeval_ortho_expr[,c(3,5)],1,FUN=function(x){max(x,na.rm=T)})
}else
{
  if  (combmethod=="mean")
  {
    sub_finalresforeval_ortho_expr$infWcom<-apply(sub_finalresforeval_ortho_expr[,c(3,5)],1,FUN=function(x){mean(x,na.rm=T)})
  }
}


labelinf=sub_finalresforeval_ortho_expr$infWcom
truelabel=sub_finalresforeval_ortho_expr$W

AUPRCcomb=MLmetrics::PRAUC(labelinf,truelabel)
print(paste0("The AUPR score for expression + ortho is ",AUPRCcomb," we use combination with ",combmethod))


PRAUROCcomb=MLmetrics::AUC(labelinf,truelabel)
print(paste0("The AUROC score for expression + ortho is ",PRAUROCcomb," we use combination with ",combmethod))
pred <- prediction(labelinf,truelabel)
roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")

plot(roc.perf,main="ROCR for expression data + orthology",colorize=T)
abline(a=0, b= 1)
auc.perf = performance(pred, measure = "auc")
auc.perf@y.values
pr.perf = performance(pred, measure="prec", x.measure="rec")
plot(pr.perf,main="PRCR for expression data + prior",colorize=T)
score=modcalcAUPRC(performance = pr.perf, subdivisions = 10000,method = "integrate")
print(score)
x<-pr.perf@x.values[[1]]
x<-x[2:length(x)]
y<-pr.perf@y.values[[1]]
y<-y[2:length(y)]
pracma::trapz(x,y)



pROC_obj <-pROC::roc(truelabel,labelinf,smoothed = TRUE,
                     # arguments for ci
                     ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                     # arguments for plot
                     plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                     print.auc=TRUE, show.thres=TRUE)

#labelinf=sub_finalresforeval_ortho_expr$infWcom
#truelabel=sub_finalresforeval_ortho_expr$W
fg <- subset(sub_finalresforeval_ortho_expr,W==1)$infWcom
bg <- subset(sub_finalresforeval_ortho_expr,W==0)$infWcom
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
plot(pr)

roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
plot(roc)

# # ############################ Rattus
# 
# orthologdata_from_rattus<-subset(orthologdata,Species%in%"10116")
# #genename<-unlist(lapply(orthologdata_from_zebra$Query,FUN=function(x){strsplit(x,"|",fixed = T)[[1]][3]}))
# #genename<-unlist(lapply(genename,FUN=function(x){strsplit(x,"_",fixed = T)[[1]][1]}))
# #orthologdata_from_zebra$genename<-genename
# multannot<-apply(orthologdata_from_rattus,1,FUN=selectrowmultannot)
# multannot<-multannot[!unlist(lapply(multannot, is.null))]
# 
# for (annotinfo in multannot)
# {
#   nbrow<-nrow(orthologdata_from_rattus)
#   pos=nbrow+1
#   listannot<-unlist(strsplit(annotinfo[4],",",fixed = T)[[1]])
#   nbannot<-length(listannot)
#   
#   if (nbannot>1)
#   {
#     uniprotID<-annotinfo[1]
#     query<-annotinfo[2]
#     species<-annotinfo[3]
#     genename<-annotinfo[5]
#     for (indx in 1:nbannot)
#     {
#       ortholog<-listannot[indx]
#       #print(ortholog)
#       orthologdata_from_rattus[pos,]<-c( uniprotID,query,species,ortholog,genename)
#       #pos=pos+1
#     }
#   }
# }
# orthologdata_from_rattus<-orthologdata_from_rattus[!str_detect(orthologdata_from_rattus$Orthologs,","),]
# orthologdata_from_rattus$ID.ortholog<-unlist(lapply(orthologdata_from_rattus$Orthologs,FUN=function(x){strsplit(x,".",fixed = T)[[1]][2]}))
# 
# 
# ############### Regulatory interaction for rattus
# fileratregnetwork<-"../data/rat_info_regulation/10116.protein.actions.v11.0.txt"
# ratregnetwork<-fread(fileratregnetwork,header = T, stringsAsFactors = F,sep="\t",
#                         colClasses = c("character","character","character","character","character","character","numeric"))
# ratregnetwork<-ratregnetwork%>%subset(mode%in%c("binding","expression") & score>=200 &item_id_a%in%orthologdata_from_rattus$Orthologs&
#                                               item_id_b%in%orthologdata_from_rattus$Orthologs&(is_directional%in%"t"|a_is_acting%in%"t"))
# 
# sub_ratregnetwork<-ratregnetwork[,c("item_id_a","item_id_b")]
# colnames(sub_ratregnetwork)<-c("TF","TG")
# sub_ratregnetwork$TF<-unlist(lapply(sub_ratregnetwork$TF,FUN=function(x){strsplit(x,".",fixed = T)[[1]][2]}))
# sub_ratregnetwork$TG<-unlist(lapply(sub_ratregnetwork$TG,FUN=function(x){strsplit(x,".",fixed = T)[[1]][2]}))
# 
# ### Potential tf
# rattffile<-"../data/rat_info_regulation/Rattus_norvegicus_TF.txt"
# rattfinfo<-read.table(rattffile,stringsAsFactors = F,header = T, quote="",sep="\t",fill=T)
# sub_rattfinfo<-rattfinfo[,c("Symbol","Ensembl", "Family", "Protein", "Entrez.ID")]
# multannot<-apply(sub_rattfinfo,1,FUN=selectrowmultannot,sep=";")
# multannot<-multannot[!unlist(lapply(multannot, is.null))]
# genewithmultprot<-unlist(lapply(multannot,FUN=function(x)x[[1]]))
# clean_sub_rattfinfo<-subset(sub_rattfinfo,!Symbol%in%genewithmultprot)
# # 
# # 
# for (annotinfo in multannot)
# {
#   nbrow<-nrow(clean_sub_rattfinfo)
#   pos=nbrow+1
#   listannot<-unlist(strsplit(annotinfo[4],";",fixed = T)[[1]])
#   nbannot<-length(listannot)
#   if (nbannot>1)
#   {
# 
#     Symbol<-annotinfo[1]
#       Ensembl<-annotinfo[2]
#       Family<-annotinfo[3]
#     Entrez.ID<-annotinfo[5]
#     for (indx in 1:nbannot)
#     {
#       protein<-listannot[indx]
#       if (protein!="")
#       {
#         #print(c(Species,Symbol,protein,Entrez.ID))
#         clean_sub_rattfinfo[pos,]<-c(Symbol,Ensembl,Family,protein,Entrez.ID)
#       }
# 
# 
#       #print(ortholog)
# 
#       #pos=pos+1
#     }
#   }
# }
# clean_sub_rattfinfo$Protein<-gsub(";","",clean_sub_rattfinfo$Protein)
# sub_clean_sub_rattfinfo<-subset(clean_sub_rattfinfo,Protein%in%orthologdata_from_rattus$ID.ortholog)
# # 
# #sub_ratregnetwork<-subset(sub_ratregnetwork,TF%in%sub_clean_sub_rattfinfo$Protein)
# #  
# 
# 
# ortho_human_network_from_rattus<-merge(sub_ratregnetwork,orthologdata_from_rattus,
#                                       by.x="TF",by.y="ID.ortholog", all.x=T)
# final_ortho_human_network_from_rattus<-merge(ortho_human_network_from_rattus,orthologdata_from_rattus,
#                                             by.x="TG",by.y="ID.ortholog", all.x=T,suffixes=c(".TF",".TG"))
# 
# final_ortho_human_network_from_rattus<-final_ortho_human_network_from_rattus%>%subset(!is.na(genename.TF)&!is.na(genename.TG))
# 
# final_ortho_human_network_from_rattus<-final_ortho_human_network_from_rattus[,c("genename.TF","genename.TG")]
# 
# colnames(final_ortho_human_network_from_rattus)<-c("TF","TG")
# 
# #### HMG1 is HMGB1
# final_ortho_human_network_from_rattus[final_ortho_human_network_from_rattus$TF%in%"HMG1","TF"]<-"HMGB1"
# #test<-sub_zebranetwork%>%subset(TF%in%sub_TFInfomotifzebra$Protein)
# sub_final_ortho_human_network_from_rattus<-subset(final_ortho_human_network_from_rattus,TF%in%intersect(humantfinfo$Name,cell_cycles_genes))
# sub_final_ortho_human_network_from_rattus<-sub_final_ortho_human_network_from_rattus%>%distinct(TF,TG,.keep_all =T)
# 
# sub_final_ortho_human_network_from_rattus$infW<-1
# 
# info_network_ortho_expr_rattus<-merge(allfinalinferrednetwork[,c("TF","TG","infW")],
#                                      sub_final_ortho_human_network_from_rattus,by=c("TF","TG"),sort=F,all.y=T,suffixes=c(".expr",".orth"))
# info_network_ortho_expr_rattus<-merge(uniquegold_standard_network,info_network_ortho_expr_rattus,by=c("TF","TG"),sort=F,all.y=T)
# 
# color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
# ######choice from 334, 26,142,74,232,329,411
# 
# 
# ###### building network
# #network_gs<-graph_from_data_frame(info_network_ortho_expr[!is.na(W)&W>0,c("TF","TG","W")])
# network_gs<-graph_from_data_frame(subset(info_network_ortho_expr_rattus,!is.na(W)&W>0,select=c("TF","TG","W")),directed = T)
# #darkcyan
# network_gs<-network_gs%>%set_edge_attr("color", value = color[74])
# #graph_network_expr<-graph_from_data_frame(info_network_ortho_expr[!is.na(infW.expr)& infW.expr>=0.5,c("TF","TG","infW.expr")],directed = T)
# graph_network_expr<-graph_from_data_frame(subset(info_network_ortho_expr_rattus,!is.na(infW.expr)& infW.expr>=0.5
#                                                  ,select=c("TF","TG","infW.expr")),directed = T)
# 
# 
# #purple
# graph_network_expr<-graph_network_expr%>%set_edge_attr("color", value = color[329])
# #graph_network_ortho<-graph_from_data_frame(info_network_ortho_expr[!is.na(infW.orth),c("TF","TG","infW.orth")],directed = T)
# graph_network_ortho<-graph_from_data_frame(subset(info_network_ortho_expr_rattus,!is.na(infW.orth),select=c("TF","TG","infW.orth")),directed = T)
# 
# #green
# graph_network_ortho<-graph_network_ortho%>%set_edge_attr("color", value = color[142])
# allnetwork<-graph_network_ortho+graph_network_expr+network_gs
# nbedges<-length(E(allnetwork)$color)
# edgecol<-unlist(lapply(1:nbedges, FUN =getcolor,E(allnetwork)$color,E(allnetwork)$color_1,E(allnetwork)$color_2,color ))
# E(allnetwork)$color<-edgecol
# ### Saving the orthonetwork from mouse to a file
# model="rattus"
# namefilegraph=paste0("../res_in_sillico/res_orthology/res_from_orthology_with_",model,"_withexpr.gml")
# #write_graph(allnetwork,namefilegraph,"gml")
# 
# 
# ################ Combining
# allinffedgesorth<-merge(info_network_ortho_expr_rattus,info_network_ortho_expr_zebra,by=c("TF","TG"),all=T,suffixes=c(".rattus",".zebra"))
# allinffedgesorth<-merge(allinffedgesorth,info_network_ortho_expr,by=c("TF","TG"),all=T)
# allinffedgesorth<-allinffedgesorth[,c("TF", "TG", "infW.orth.rattus","infW.orth.zebra","infW.orth","W", "infW.expr")]
# sub_allinffedgesorth<-subset(allinffedgesorth,infW.expr>=0.5)
# intersection<-subset(sub_allinffedgesorth,!is.na(infW.orth.rattus)&!is.na(infW.orth.zebra)&!is.na(infW.orth)&!is.na(W))
# intersection_ratzebra<-subset(sub_allinffedgesorth,!is.na(infW.orth.rattus)&!is.na(infW.orth.zebra)&is.na(infW.orth)&is.na(W))
# intersection_ratmouse<-subset(sub_allinffedgesorth,!is.na(infW.orth.rattus)&is.na(infW.orth.zebra)&!is.na(infW.orth)&is.na(W))
# intersection_zebramouse<-subset(sub_allinffedgesorth,is.na(infW.orth.rattus)&!is.na(infW.orth.zebra)&!is.na(infW.orth)&is.na(W))
# intersection_new<-subset(sub_allinffedgesorth,!is.na(infW.orth.rattus)&!is.na(infW.orth.zebra)&!is.na(infW.orth)&is.na(W))

######### Ploting venndidgram
#myV <- createVennObj(nSets = 5, sNames = c('expr', 'orth', 'gs-NA', 'gs-True', 'gs-False'), 
 #                    sSizes = c(0,156,44,0,343,0,0,0,543,156,44,0,343,0,0,0,403,143,39,0,221,0,0,0,403,145,39,0,221,0,0,0))


#myvennimg<-plotVenn(nVennObj = myV, setColors = c("skyblue", "pink1", "mediumorchid", "yellow", "orange"), borderWidth = 0,
#                    outFile="../res_in_sillico/res_orthology/venn-ortho-plot.svg",systemShow =T)







# 
# fileratregnetwork2<-"../data/rat_info_regulation/R_norvegicus_interactions.txt"
# ratregnetwork2<-fread(fileratregnetwork2,header = F, stringsAsFactors = F,sep="\t",
#                       colClasses = c("character","character","character","character","character","character",
#                                      "character","character","character","character","character","numeric","numeric","numeric","numeric","character"),skip=5)
# 


####################Trash



# R_norvegicus_interactions <- read_delim("/Volumes/Seagate Backup Plus Drive/PhD/GRN_inference_prime/data/rat_info_regulation/R_norvegicus_interactions.txt", 
#                                       "\t", escape_double = FALSE, trim_ws = TRUE)
# 
# fileratregnetwork2<-"../data/rat_info_regulation/interaction_rat.txt"
# ratregnetwork2<-fread(fileratregnetwork2,header = F, stringsAsFactors = F,sep="\t",colClasses = c("character","character","character","character","character",
#                                                                                                 "character","character","character","character","character","character"),skip=10)
# ratregnetwork2<-ratregnetwork2[,c(1:10)]
# colnames(ratregnetwork2)<-c("Interactor_A", 	"Interactor_A_Gene",	"Interactor_A_Gene RGD_ID","Interactor_B","Interactor_B_Gene",
#                            "Interactor_B_Gene_RGD_ID",	"Species_A",	"Species_B",	"Interaction_Type",	"Attributes")
# sub_ratregnetwork2<-subset(ratregnetwork2,ratregnetwork2$Species_A=="Rat"&ratregnetwork2$Species_B=="Rat"&ratregnetwork2$Interaction_Type=="physical association")
# # filemappinguniprotid<-"../data/rat_info_regulation/mapping_prot_string.txt"
# mappinguniprotid<-read.table(filemappinguniprotid,stringsAsFactors = F,header = T,sep="\t")
# mappinguniprotid$Protein<-unlist(lapply(mappinguniprotid$To,FUN=function(x){strsplit(x,".",fixed = T)[[1]][2]}))
# final_sub_ratregnetwork2<-merge(sub_ratregnetwork2,mappinguniprotid,by.x="Interactor_A",by.y="From")
# final_sub_ratregnetwork2<-merge(final_sub_ratregnetwork2,mappinguniprotid,by.x="Interactor_B",by.y="From",suffixes=c(".TF",".TG"))
# final_sub_ratregnetwork2<-final_sub_ratregnetwork2[,c("Protein.TF","Protein.TG")]
# test<-subset(final_sub_ratregnetwork2,Protein.TF%in%orthologdata_from_rattus$ID.ortholog&Protein.TG%in%orthologdata_from_rattus$ID.ortholog)
# 
# 
# 
# 

# 
# ####### 
# filenetworkrat<-""
# 
# 
# 
# final_mouseregulatorynetwork<-rbind(sub_mouseregnetwork,sub_mouseregnetwork2)
# final_mouseregulatorynetwork<-
# unique_mouseregulatorynetwork<-final_mouseregulatorynetwork%>%distinct(Gene1,Gene2,.keep_all =T)
# 
# 
# 
# 
# 
# #BENINKozscore(network=1,alphaenet=0.9,nbBoobstrap=1000,exponent=1.4,sizenetwork=10,lambda=20,lmean=10,alpha=0.5,lambdamin=1,lambdamax=1000,potentialtf="known")
# # 
# # source("bootfuns.q")
# # source("bootpracs.q")
# # library("BiocGenerics",lib.loc="./library")
# # library("backports",lib.loc="./library")
# # library("pkgconfig",lib.loc="./library")
# # 
# # BENINKozscore(network=5,alphaenet=0.9,nbBoobstrap=4000,exponent=1.4,sizenetwork=100,lambda=20,lmean=10,alpha=0.5,lambdamin=1,lambdamax=1000,potentialtf="all")


mouseortho<-igraph::read_graph("/Users/Stephanie/Downloads/res_from_orthology_with_mouse_new.gml",format="gml")
mydataframe<-as_data_frame(mouseortho, what = c("both"))
edgesmouseortho<-mydataframe$edges
foxm1tgorth<-subset(edgesmouseortho,from%in%"FOXM1")




#clustering
exprdatfile<-"../data/data_human/final_data_hum_reg_network/Hela_data/imputed_data_human_unique_ens_id.txt"

final_annot_human_exprdata_unique_id<-read.table(exprdatfile,header=T, stringsAsFactors =F, fill=T,quote="",sep="\t")
row.names(final_annot_human_exprdata_unique_id)<-final_annot_human_exprdata_unique_id$Symbol
final_annot_human_exprdata_unique_id<-final_annot_human_exprdata_unique_id[,-1]
considered_exprdat<-final_annot_human_exprdata_unique_id[,39:86]
#write.table(row.names(considered_exprdat),"../data/data_human/final_data_hum_reg_network/Hela_data/considered_genes_hela_cell_cycle.txt",col.names = F,
#            row.names = F,quote=F)
scaled_considered_exprdat<-scale(considered_exprdat)
corr_hela_expr_data<- cor(t(scaled_considered_exprdat))
gene_dist<-dist(corr_hela_expr_data)
gene_hclust <- hclust(gene_dist)
mycol <- redgreen(75)
mycl <- cutree(gene_hclust, k = 5)
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 

heatmap.2(scaled_considered_exprdat, Rowv=as.dendrogram(hr), col=mycol, scale="row", density.info="none", trace="none", RowSideColors=mycolhc)

