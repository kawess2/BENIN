setwd("/Users/stephanie_admin/Documents/GRN_inference/src/")
library(tidyverse)
library(dplyr)
library(reshape2)
gsnetwork<-read.csv("../data/dream4_challenge_data/gold_standard/DREAM4_Challenge2_GoldStandards/DREAM4_GoldStandard_InSilico_Size100_5.tsv",sep = '\t',stringsAsFactors =F)
colnames(gsnetwork)<-c("TF","TG","W")
network_gs<-graph_from_data_frame(gsnetwork,directed = T)
subgsnetwork<-subset(gsnetwork,TF=="G74")
E(subgsnetwork)$color="black"
g1=induced_subgraph(network_gs,vids=unique(subgsnetwork$TF,subgsnetwork$TG),impl ="copy_and_delete")
networkbenin<-read.csv("../res_in_sillico/insilico_size100_5/lambda_20/sub_network_network_5_benin.tsv",sep = '\t',stringsAsFactors =F)
colnames(networkbenin)<-c("TF","TG","W")
gs_diff<-full_join(gsnetwork,networkbenin,by=c("TF","TG"),suffix=c(".gs",".benin"))
gs_diff[is.na(gs_diff)]<-0
network_benin<-graph_from_data_frame(networkbenin,directed = T)
network_benin<-network_benin%>%set_edge_attr("color", value = "green")
g2<-induced_subgraph(network_benin,vids=unique(union(subnetworkbenin$TF,subnetworkbenin$TG)),impl ="copy_and_delete")
g<-g1+g2
E(g)$color=ifelse(is.na(E(g)$color_2), E(g)$color_1,E(g)$color_2)
g3<-g2-g1
E(g3)$color="green"
gprime=g+g3
E(gprime)$color=ifelse(is.na(E(gprime)$color_2), E(gprime)$color_1,E(gprime)$color_2)
plot(gprime)
# listtf<-unique(as.character(networkbenin$TF))
# listtg<-unique(as.character(networkbenin$TG))
# networkbenin<-networkbenin%>%dcast(TG~TF,fill=0)
# row.names(networkbenin)<-networkbenin[,1]
# networkbenin<-networkbenin[,-1]
# networkbenin<-networkbenin[listtg,listtf]
networkscanbma<-read.csv("../res_in_sillico/insilico_size100_5/lambda_20/sub_network_network_5_scanBMA.tsv",sep = '\t',stringsAsFactors =F)
colnames(networkscanbma)<-c("TF","TG","W")
gs_diff<-full_join(gs_diff,networkscanbma,by=c("TF","TG"))
gs_diff[is.na(gs_diff)]<-0
networkirafnet<-read.csv("../res_in_sillico/insilico_size100_5/lambda_20/sub_network_100_5_irafnet.tsv",sep = '\t',stringsAsFactors =F)
colnames(networkirafnet)<-c("TF","TG","W")
gs_diff<-full_join(gs_diff,networkirafnet,by=c("TF","TG"),suffix=c(".scanbma",".irafnet"))
gs_diff[is.na(gs_diff)]<-0
networkgelnet<-read.csv("../res_in_sillico/insilico_size100_5/lambda_20/sub_network_100_5_gelnet.tsv",sep = '\t',stringsAsFactors =F)
colnames(networkgelnet)<-c("TF","TG","W")
gs_diff<-full_join(gs_diff,networkgelnet,by=c("TF","TG"),suffix=c(".scanbma",".irafnet"))
gs_diff[is.na(gs_diff)]<-0
nbrow=dim(gs_diff)[1]
for(i in 1:nbrow)
{
  print(isTRUE(gs_diff[i,c("W.gs","W.benin","W.scanbma","W.irafnet","W")]))
}
# > TU_col =edge_attr(g, "color_2")
# > E2 = which(is.na(edge_attr(tree_union, "color_1")))
# Error in "igraph" %in% class(graph) : object 'tree_union' not found
# > E2 = which(is.na(edge_attr(g, "color_1")))
# > E2 = which(is.na(edge_attr(g, "color_2")))
# > TU_col[E2] = edge_attr(g, "color_2")[E2]
# > TU_col[E2] = edge_attr(g, "color_1")[E2]
# > tree_union2 = set_edge_attr(g, "color", value=TU_col)
# > plot(tree_union2)
