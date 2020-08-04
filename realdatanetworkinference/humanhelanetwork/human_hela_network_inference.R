#setwd("/Users/stephanie_admin/Documents/GRN_inference/src/")
library(glmnet)
library(xts)
library(boot)
library(caTools)
library(RGeode)
library(dplyr)
library(doSNOW)
library(parallel)
library(foreach)
library(snow)
library(iterators)
library(zoo)
library(xts)
library(boot)
library(caTools)
library(crayon)
library(data.table)
library(reshape2)
library(doRNG)
#library(TFTargetCaller)
#library(universalmotif)
library(DMwR)
library(conclust)
library(tidyverse)
library(data.table)
library(seqinr)
library(STRINGdb)
library(GOSemSim)
library(org.Hs.eg.db)
library(readr)

source("R/bootfuns.q")
source("R/bootpracs.q")
source("R/benin.R")
source("R/utile.R")




####################################
#We first analyze the expression data for clustering to get the module network

# human_expr_file<-"../data/data_human/final_data_hum_reg_network/Hela_data/hela_ts_data.txt"
# all_human_expr_data<-read.table(human_expr_file,header=T, stringsAsFactors =F, fill=T,quote="",sep="\t")
# all_human_expr_data<-all_human_expr_data[, colSums(is.na(all_human_expr_data)) != nrow(all_human_expr_data)]
# all_human_expr_data<-all_human_expr_data[,c(1,13:126)]
# row.names(all_human_expr_data)<-all_human_expr_data$UID
# all_human_expr_data<-all_human_expr_data[,-1]
# gene_cell_cycle_hela_file<-"../data/data_human/final_data_hum_reg_network/Hela_data/cell-cycle-hela-data-whitfield.txt"
# gene_cell_cycle_hela<-read.table(gene_cell_cycle_hela_file,header=T, stringsAsFactors =F, fill=T,quote="",sep="\t")
# sub_human_expr_data<-all_human_expr_data[gene_cell_cycle_hela$CloneID,]
# sub_gene_cell_cycle_hela<-subset(gene_cell_cycle_hela,!(gene_cell_cycle_hela$Symbol==""))
# sub_human_expr_data<-all_human_expr_data[gene_cell_cycle_hela$CloneID,]
# final_annot_human_exprdata_mult_id<-knnImputation(sub_human_expr_data,k=12)
# merged_final_annot_human_exprdata_mult_id<-merge(final_annot_human_exprdata_mult_id,gene_cell_cycle_hela[,c(1,4)],by.x=0,by.y="CloneID")
# merged_final_annot_human_exprdata_unique_id<-data.table(merged_final_annot_human_exprdata_mult_id)
# final_annot_human_exprdata_unique_id<-(subset(merged_final_annot_human_exprdata_mult_id,!(merged_final_annot_human_exprdata_mult_id$Symbol=="")))
# final_annot_human_exprdata_unique_id<-as.data.table(final_annot_human_exprdata_unique_id)
# final_annot_human_exprdata_unique_id<-final_annot_human_exprdata_unique_id[,-1]
# 
# final_annot_human_exprdata_unique_id<-final_annot_human_exprdata_unique_id[, lapply(.SD,mean),by = Symbol]
# write.table(final_annot_human_exprdata_unique_id, file="../data/data_human/final_data_hum_reg_network/Hela_data/imputed_data_human_unique_ens_id.txt",row.names = F,quote=F,sep= "\t")
# 
# #clustering
# exprdatfile<-"../data/data_human/final_data_hum_reg_network/Hela_data/imputed_data_human_unique_ens_id.txt"
# 
# final_annot_human_exprdata_unique_id<-read.table(exprdatfile,header=T, stringsAsFactors =F, fill=T,quote="",sep="\t")
# row.names(final_annot_human_exprdata_unique_id)<-final_annot_human_exprdata_unique_id$Symbol
# final_annot_human_exprdata_unique_id<-final_annot_human_exprdata_unique_id[,-1]
# considered_exprdat<-final_annot_human_exprdata_unique_id[,39:86]
# #write.table(row.names(considered_exprdat),"../data/data_human/final_data_hum_reg_network/Hela_data/considered_genes_hela_cell_cycle.txt",col.names = F,
# #            row.names = F,quote=F)
# scaled_considered_exprdat<-scale(considered_exprdat)
# corr_hela_expr_data<- cor(t(scaled_considered_exprdat))
# gene_dist<-dist(corr_hela_expr_data)
# gene_hclust <- hclust(gene_dist)
# plot(gene_hclust, labels = FALSE)
# gene_cluster <- cutree(gene_hclust, k = 5) %>% 
#   # turn the named vector into a tibble
#   enframe() %>% 
#   # rename some of the columns
#   rename(gene = name, cluster = value)
# write.table(subset(gene_cluster,cluster==1)$gene, file="../data/data_human/final_data_hum_reg_network/Hela_data/modules/cluster1_comp_clust.txt",row.names = F,
#             quote=F,sep= "\t",col.names = F)
# write.table(subset(gene_cluster,cluster==2)$gene, file="../data/data_human/final_data_hum_reg_network/Hela_data/modules/cluster2_comp_clust.txt",row.names = F,
#             quote=F,sep= "\t",col.names = F)
# write.table(subset(gene_cluster,cluster==3)$gene, file="../data/data_human/final_data_hum_reg_network/Hela_data/modules/cluster3_comp_clust.txt",row.names = F,
#             quote=F,sep= "\t",col.names = F)
# write.table(subset(gene_cluster,cluster==4)$gene, file="../data/data_human/final_data_hum_reg_network/Hela_data/modules/cluster4_comp_clust.txt",row.names = F,
#             quote=F,sep= "\t",col.names = F)
# write.table(subset(gene_cluster,cluster==5)$gene, file="../data/data_human/final_data_hum_reg_network/Hela_data/modules/cluster5_comp_clust.txt.txt",row.names = F,
#             quote=F,sep= "\t",col.names = F)
# write.table(subset(gene_cluster,cluster==6)$gene, file="../data/data_human/final_data_hum_reg_network/Hela_data/modules/cluster6_comp_clust.txt",row.names = F,
#             quote=F,sep= "\t",col.names = F)
# 
# #Get the ens_transcript ids for the genes in clusters to get their promoter sequences
# mapping_genes_from_uscs<-read.table("../data/data_human/final_data_hum_reg_network/Hela_data/sequences/mapping_gene_uscs_genome.txt",header = F,
#                                     stringsAsFactors = F,fill = T,sep="\t")
# sub_mapping_genes_from_uscs<-subset(mapping_genes_from_uscs,mapping_genes_from_uscs$V3%in%row.names(final_annot_human_exprdata_unique_id))
# 
# write.table(sub_mapping_genes_from_uscs$V1,"../data/data_human/final_data_hum_reg_network/Hela_data/ens_id.txt",
#             sep = "\t",row.names = F,col.names = F,quote = F)
# # We split the promoter sequences in different files corresponding to each clusters
# mypromoter<-read.fasta(file = "../data/data_human/final_data_hum_reg_network/Hela_data/sequences/promoter_sequences_all.fa",
#                        seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
# 
# #hg38_knownGene_
# nbcluster=1
# for (id_cluster in 1:nbcluster)
# {
#   name_file_seq_in_cluster<-paste0("../data/data_human/final_data_hum_reg_network/Hela_data/sequences/cluster_",id_cluster ,".fa")
#   #id_genes_cluster<-unlist(lapply(subset(sub_mapping_genes_from_uscs,sub_mapping_genes_from_uscs$V3%in%subset(gene_cluster,cluster==id_cluster)$gene)$V1,
#   #              FUN=function(x,y)paste0(y,x),y="hg38_knownGene_"))
#   id_genes_cluster<-unlist(lapply(subset(sub_mapping_genes_from_uscs,sub_mapping_genes_from_uscs$V3%in%subset(gene_cluster,cluster==id_cluster)$gene)$V1,
#                 FUN=function(x,y)paste0(y,x),y=""))
#   sequences_in_files<-mypromoter[id_genes_cluster]
#   namesseq<-names(sequences_in_files)
#   write.fasta(sequences_in_files,namesseq,name_file_seq_in_cluster)
# }
# 
# ## Getting the list of transcription factors
# fileTFInfomotif<-paste0("/Users/stephanie_admin/Documents/GRN_inference/data/data_human/final_data_hum_reg_network/",
#                    "Hela_data/Homo_sapiens_2020_02_24_4-34_pm/TF_Information.txt")
# TFInfomotif<-read.table(fileTFInfomotif, header = T, stringsAsFactors = F, sep="\t",fill = T)
# filepotentialTF<-paste0("/Users/stephanie_admin/Documents/GRN_inference/data/data_human/final_data_hum_reg_network/Hela_data/human_TF_list.txt")
# potentialTF<-read.table(filepotentialTF, header = T, stringsAsFactors = F, sep="\t",fill = T)
# potentialTF<-subset(potentialTF,potentialTF$Is_TF=="Yes")
# potentialTFname<-intersect(sub_gene_cell_cycle_hela$Symbol,TFInfomotif$TF_Name)
# write.table(potentialTFname,"../data/data_human/final_data_hum_reg_network/Hela_data/Hela_cell_cycle_TF.txt",col.names = F, row.names = F,quote=F)
# # 2) First we get the list of motifs file
# mymotifsfilelist=list.files("/Users/stephanie_admin/Documents/GRN_inference/data/data_human/final_data_hum_reg_network/Hela_data/Homo_sapiens_2020_02_24_4-34_pm/pwms_all_motifs",no.. = T)
# listmotifid=lapply(mymotifsfilelist,FUN=function(x){gsub('.txt','',x)})
# #3 ) Get the the description of of the list of motif in the directory and which are from Transfac, chip-seq, chip-chip
# #sub_description_motif=subset(description_motif,(Motif_ID%in%listmotifid & Motif_Type %in%c("Misc","Transfac","ChIP-seq","ChIP-chip")))
# sub_description_motif<-subset(TFInfomotif,Motif_ID !=".")
# sub_description_motif<-subset(sub_description_motif,sub_description_motif$TF_Name %in%potentialTFname)
# #We finaly get #numbers of TF
# nbtfwmotif<-length(unique (sub_description_motif$TF_Name))
# # For a total of #motif
# nbmotif<-dim(sub_description_motif)[1]
# #We need to transform the motif in meme format
# allmotifsfilename="/Users/stephanie_admin/Documents/GRN_inference/data/data_human/final_data_hum_reg_network/Hela_data/Homo_sapiens_2020_02_24_4-34_pm/allmofifscisbp.txt"
# lapply(seq(1,nbmotif),writemotif, allmotifsfilename=allmotifsfilename,matallmotif=sub_description_motif)
# cisbpmotifs<-read_cisbp(allmotifsfilename)
# memecisbpmotifsfilename="../data/data_human/final_data_hum_reg_network/Hela_data/Homo_sapiens_2020_02_24_4-34_pm/Homo_sapiens.meme"
# write_meme(cisbpmotifs,memecisbpmotifsfilename)
# 
# 
# # Go semantic similarity
# hsGO<-godata("org.Hs.eg.db",keytype="SYMBOL",ont = "BP",computeIC = TRUE)
# genes<-unique(sub_gene_cell_cycle_hela$Symbol)
# semanticsim <-mgeneSim(genes,semData = hsGO,measure ="Rel",combine="BMA" )
# copy_semanticsim<-semanticsim
# #hsGO2<-godata("org.Hs.eg.db",keytype="ENTREZID",ont = "BP",computeIC = TRUE)
# #semanticsim2<-mgeneSim(unique(sub_gene_cell_cycle_hela$EntrezID_W),semData = hsGO2,measure ="Rel",combine="BMA" )
# #******************************************************************************************************************************
# 
# #We now built features for eachcluster to do logistic regression
# 
# # hela_gs_network<-full_join(hela_gs_network1,hela_gs_network2,by=c("From","To"))
# # final_hela_gs_network<-left_join(hela_gs_network,mapping_tf_entrez_name,by=c("From"="ENTREZID"))
# # final_hela_gs_network<-left_join(final_hela_gs_network,mapping_entrez_name,by=c("To"="ENTREZID"))
# # final_hela_gs_network<-final_hela_gs_network[,c("SYMBOL.x","SYMBOL.y","W")]
# # colnames(final_hela_gs_network)<-c("TF","TG","W")
# #hela_gs_network[hela_gs_network$W!=1,]<-0
# #hela_gs_network<-hela_gs_network[order(hela_gs_network$W,decreasing = T),]
# #hela_gs_network$From<-as.character(hela_gs_network$From)
# #hela_gs_network$To<-as.character(hela_gs_network$To)
# #mapping_entrez_name<-select(org.Hs.eg.db,as.character(unique(c(hela_gs_network$From,hela_gs_network$To))), c("ENTREZID","GENENAME","ENSEMBL","SYMBOL"), "ENTREZID")
# ################################ Code for building our human gold standar
# ##Getting hela cell cycle genes
# file_cell_cycle_genes<-"../data/data_human/final_data_hum_reg_network/Hela_data/considered_genes_hela_cell_cycle.txt"
# cell_cycles_genes<-(read.table(file_cell_cycle_genes,header = F,quote="",sep="\t",stringsAsFactors = F))$V1
# file_potential_tf<-"../data/data_human/final_data_hum_reg_network/Hela_data/Hela_cell_cycle_TF.txt"
# potentialTFname<-as.vector(read.table(file_potential_tf,stringsAsFactors = F,quote = "",header = F)$V1)
# mappingnamefile<-"../data/data_human/final_data_hum_reg_network/Hela_data/cell-cycle-hela-data-whitfield.txt"
# mapping_entrez_name2<-read.table(mappingnamefile,header = T,quote="",sep="\t",stringsAsFactors = F,fill=T,colClasses = c("character","character","character","character","character","character"))
# mapping_entrez_name2<-subset(mapping_entrez_name2,mapping_entrez_name2$Symbol%in%cell_cycles_genes)
# mapping_entrez_name2<-mapping_entrez_name2[,c("EntrezID_W","Symbol")]
# colnames(mapping_entrez_name2)<-c("ENTREZID","SYMBOL")
# mapping_tf_entrez_name2<-subset(mapping_entrez_name2,mapping_entrez_name2$SYMBOL%in%potentialTFname)
# #### Mapping from org.Hs.eg.db
# mapping_entrez_name<-select(org.Hs.eg.db,unique(cell_cycles_genes),
#                             c("ENTREZID","GENENAME","ENSEMBL","SYMBOL"), "SYMBOL")
# mapping_tf_entrez_name<-select(org.Hs.eg.db,unique(potentialTFname),
#                                c("ENTREZID","GENENAME","ENSEMBL","SYMBOL"), "SYMBOL")
# 
# notmappedgenes<-inner_join(mapping_entrez_name[is.na(mapping_entrez_name$ENTREZID),],mapping_entrez_name2,by="SYMBOL")
# notmappedgenes<-unique(notmappedgenes)
# mapping_entrez_name[mapping_entrez_name$SYMBOL%in%notmappedgenes$SYMBOL,"ENTREZID"]<-notmappedgenes$ENTREZID.y
# #Reading golstandard network 
# #file_gs_cervix<-"../data/data_human/final_data_hum_reg_network/Hela_data/uterine_cervix_gold_standard.dat"
# #######Option 1
# #hela_gs_network<-read.table(file_gs_cervix,header = F,quote="",sep="\t",stringsAsFactors = F,col.names=c("From","To","W"))
# file_gs_hela<-"../data/data_human/final_data_hum_reg_network/Hela_data/gs_file_tissue/global_hela_gs.dat"
# hela_gs_network1<-fread(file=file_gs_hela,header = F,quote="",sep="\t",stringsAsFactors = F,col.names=c("From","To","W"),
#                         colClasses = c("character","character","numeric"))
# hela_gs_network1<-subset(hela_gs_network1,hela_gs_network1$From%in%mapping_tf_entrez_name$ENTREZID&hela_gs_network1$To%in%mapping_entrez_name$ENTREZID)
# hela_gs_network1<-left_join(hela_gs_network1,mapping_tf_entrez_name,by=c("From"="ENTREZID"))
# hela_gs_network1<-left_join(hela_gs_network1,mapping_entrez_name,by=c("To"="ENTREZID"))
# hela_gs_network1<-hela_gs_network1[,c("SYMBOL.x","SYMBOL.y","W")]
# hela_gs_network1<-unique(hela_gs_network1)
# colnames(hela_gs_network1)<-c("TF","TG","W")
# ######Option 2
# #file_gs_hela<-"../data/data_human/final_data_hum_reg_network/Hela_data/human_gs_TFTG_network.csv"
# #Here use import it is much more easy
# #human_gs_TFTG_network
# human_gs_TFTG_network <- read_csv("~/Documents/GRN_inference/data/data_human/final_data_hum_reg_network/Hela_data/human_gs_TFTG_network.csv")
# hela_gs_network2<-subset(human_gs_TFTG_network,!human_gs_TFTG_network$is_evidence_chipSeq&!human_gs_TFTG_network$is_evidence_TFbindingMotif)
# hela_gs_network2<-subset(hela_gs_network2,hela_gs_network2$TF%in%potentialTFname &hela_gs_network2$target%in%cell_cycles_genes )
# final_human_gs<-full_join(hela_gs_network2,hela_gs_network1,by=c("TF"="TF","target"="TG"))
# final_human_gs[is.na(final_human_gs$W),"W"]<-1
# final_human_gs[(final_human_gs$W!=1),"W"]<-0
# final_human_gs<-final_human_gs[,c("TF","target","W")]
# colnames(final_human_gs)<-c("TF","TG","W")
# 
# 









################################
#Getting 
#working on cluster




networkinferencepercluster<-function(priordata,clusterid,lambdamin=1,lambdamax=1000,beta=0.5,exponent=1.1,alphaenet=0.6,nbBootstrap=1000,ont="BP",methodfuncsim="Rel",sizenetwork=0)
{
  
  # We first get the data from the cluster
  clusterid= strsplit(clusterid,"_")[[1]][2]
  file_cell_cycle_genes<-"../data/data_human/final_data_hum_reg_network/Hela_data/considered_genes_hela_cell_cycle.txt"
  cell_cycles_genes<-(read.table(file_cell_cycle_genes,header = F,quote="",sep="\t",stringsAsFactors = F))$V1
  file_potential_tf<-"../data/data_human/final_data_hum_reg_network/Hela_data/Hela_cell_cycle_TF.txt"
  potentialTFname<-as.vector(read.table(file_potential_tf,stringsAsFactors = F,quote = "",header = F)$V1)
  mappingnamefile<-"../data/data_human/final_data_hum_reg_network/Hela_data/cell-cycle-hela-data-whitfield.txt"
  mapping_entrez_name2<-read.table(mappingnamefile,header = T,quote="",sep="\t",stringsAsFactors = F,fill=T,
                                   colClasses = c("character","character","character","character",
                                                  "character","character"))
  filecluster=paste0("../data/data_human/final_data_hum_reg_network/Hela_data/modules/cluster_",clusterid,"_comp_clust.txt")
  #filecluster=paste0("../data/data_human/final_data_hum_reg_network/Hela_data/modules/",clusterid,".txt")
  genesincluster<-as.vector(read.table(filecluster,stringsAsFactors = F,header = F,quote="")$V1)
  allgenes<-unique(c(potentialTFname,genesincluster))
  numbergenes<-length(allgenes)
  
  #Reading expression data 
  exprdatfile<-"../data/data_human/final_data_hum_reg_network/Hela_data/imputed_data_human_unique_ens_id.txt"
  final_annot_human_exprdata_unique_id<-read.table(exprdatfile,header=T, stringsAsFactors =F, fill=T,quote="",sep="\t")
  row.names(final_annot_human_exprdata_unique_id)<-final_annot_human_exprdata_unique_id$Symbol
  final_annot_human_exprdata_unique_id<-final_annot_human_exprdata_unique_id[,-1]
  considered_exprdat<-final_annot_human_exprdata_unique_id[,39:86]
  filenamegs<-"../data/data_human/final_data_hum_reg_network/Hela_data/final_human_goldstandard.txt"
  human_gs<-read.table(filenamegs,header = T,stringsAsFactors = F,sep="\t")
  #Now we deal with each prior 
  #matpval<-matrixbidingprior
 
  
  matweight=matrix(0,nr = numbergenes, nc = numbergenes)
  rownames(matweight)<-allgenes
  colnames(matweight)<-allgenes
  if (priordata=="TFBS")
  {
    #analysing FIMO output
    matprior=matrix(1,nr = numbergenes, nc = numbergenes)
    
    #exponent=1.1
    matprior<-as.data.frame(matprior)
    #test<-as.data.frame(matrixbidingprior)
    rownames(matprior)<-allgenes
    colnames(matprior)<-allgenes
    resfilename=paste0("grep -v '^#' ../data/data_human/final_data_hum_reg_network/Hela_data/sequences/res_promoter_scanning/cluster_",clusterid,"/fimo.tsv")
    resfimo<-fread(resfilename, header = T, stringsAsFactors = F,sep = "\t",colClasses = c("character","character","character","character","numeric",
                                                                                           "character","numeric","numeric","numeric","character"),
                   blank.lines.skip=T)
    fileTFInfomotif<-paste0("../data/data_human/final_data_hum_reg_network/",
                            "Hela_data/Homo_sapiens_2020_02_24_4-34_pm/TF_Information.txt")
    TFInfomotif<-read.table(fileTFInfomotif, header = T, stringsAsFactors = F, sep="\t",fill = T,comment.char = "")
    sub_description_motif<-subset(TFInfomotif,Motif_ID !=".")
    sub_description_motif<-subset(sub_description_motif,sub_description_motif$TF_Name %in%potentialTFname)
    resfimo<-as.data.table(resfimo)
    sub_description_motif<-as.data.table(sub_description_motif)
    #resfimo$TGENSEMBL<-unlist(lapply(resfimo$sequence_name, FUN=function(x)strsplit(x,"_")[[1]][3]))
    resfimo$TGENSEMBL<-resfimo$sequence_name
    #setkey(resfimo, motif_id)
    #setkey(sub_description_motif, TF_ID)
    #df3 <- resfimo[sub_description_motif, nomatch = 0,allow.cartesian=TRUE]
#print("fhfjfjfjffkfkfklfklfkf")    
final_res_fimo<-inner_join(resfimo,sub_description_motif[,c(1,4,5,6,7)],by=c("motif_id"="TF_ID"))
    mapping_genes_from_uscs<-read.table("../data/data_human/final_data_hum_reg_network/Hela_data/sequences/mapping_gene_uscs_genome.txt",header = F,
                                        stringsAsFactors = F,fill = T,sep="\t")
    sub_mapping_genes_from_uscs<-subset(mapping_genes_from_uscs,mapping_genes_from_uscs$V3%in%cell_cycles_genes)
    colnames(sub_mapping_genes_from_uscs)<-c("TGENSEMBL","REFSEQTG","TG_Name","REFSEQTG2")
    #gc(resfimo)
    #gc(mapping_genes_from_uscs)
    final_res_fimo<-inner_join(final_res_fimo,sub_mapping_genes_from_uscs[,c("TGENSEMBL","TG_Name")],by="TGENSEMBL")
#print(colnames(bindingfeature))    
bindingfeature<-final_res_fimo[,c("TF_Name","DBID","TG_Name","TGENSEMBL","DBID","q-value")]
print(colnames(bindingfeature)) 
colnames(bindingfeature)<-c("TF_Name","DBID","TG_Name","TGENSEMBL","DBID","qvalue")
   bindingfeature<-bindingfeature[order(bindingfeature$qvalue,decreasing = F),]
    bindingfeature<-bindingfeature[,c("TF_Name","TG_Name","qvalue")]
#print("fhfjfjfjffkfkfklfklfkf")

    matrixbidingprior<-dcast(bindingfeature,TG_Name~TF_Name,value.var = "qvalue",fill=1,drop =F,fun.aggregate = min)
    rownames(matrixbidingprior)<-matrixbidingprior[,"TG_Name"]
    matrixbidingprior<-matrixbidingprior[,-1]
    #matpval<-matrixbidingprior
    
    #test<-as.data.frame(matrixbidingprior)
    
    matprior[rownames(matrixbidingprior),colnames(matrixbidingprior)]<-as.data.frame(matrixbidingprior)
    # Transforming location data into probabilities
    listtf<-colnames(matrixbidingprior)
    print(listtf)
    numbergenes<-length(allgenes)
    nbtf<-length(listtf)
    
    
    for (g1 in allgenes)
    {
      for (g2 in listtf)
      {
        prior=matprior[g1,g2]
        #print(pval)
        # prob=ifelse(g1==g2,0,(intergral(beta=beta,pval=pval)))
        # matweight[g1,g2]=ifelse(prob>=0.8,0.0001,1)
        matweight[g1,g2]=ifelse(g1==g2,0,1/(integral(beta=beta,pval=prior,lambdamin=lambdamin,lambdamax=lambdamax))**exponent)
      }
    }
    
    
  }
  else
  {
    if (priordata=="Chipseq")
    {
      matprior=matrix(0,nr = numbergenes, nc = numbergenes)
      
      #exponent=1.1
      matprior<-as.data.frame(matprior)
      #test<-as.data.frame(matrixbidingprior)
      rownames(matprior)<-allgenes
      colnames(matprior)<-allgenes
      ### Reading information about chip-seq data and corresponding TF
      namefilepriorchipseq<-"../data/data_human/final_data_hum_reg_network/Hela_data/chipseqdata/wgEncodeAwgTfbsUniform/TFTGscore_chipseq.txt"
      priorchipseq<-fread(namefilepriorchipseq,header = T,quote="",stringsAsFactors = F,
                        fill = T,sep="\t", colClasses = c("character","character","numeric"))
      print( unique(priorchipseq$TF))
      priorchipseq<-subset(priorchipseq,priorchipseq$TF %in%potentialTFname& priorchipseq$GeneSymbol%in%genesincluster)
      #print( unique(priorchipseq$TF))
      #print(potentialTFname)
      colnames(priorchipseq)<-c("TF","TG","W")
      matrixchipseq<-reshape2::acast(priorchipseq,TG~TF,value.var = "W",fill=0,drop =F,
                                        fun.aggregate = mean)
      
      
      matprior[rownames(matrixchipseq),colnames(matrixchipseq)]<-as.data.frame(matrixchipseq)
      # Transforming functional similarity into weight
      listtf<-colnames(matrixchipseq)
      for (g1 in allgenes)
      {
        for (g2 in listtf)
        {
          prior<- matprior[g1,g2]
          if (prior!=0) 
          {	#print (paste0(g1,"-->",g2))
            matweight[g1,g2]<-ifelse(g1==g2,0,1/(prior**exponent))
          }
        }
      }
      
    }
    else
      if (priordata=="functional")
      {
        matprior=matrix(0,nr = numbergenes, nc = numbergenes)
        
        #exponent=1.1
        matprior<-as.data.frame(matprior)
        #test<-as.data.frame(matrixbidingprior)
        rownames(matprior)<-allgenes
        colnames(matprior)<-allgenes
        hsGO<-godata("org.Hs.eg.db",keytype="SYMBOL",ont = ont,computeIC = TRUE)
        genes<-unique(allgenes)
        semanticsim_in_cluster<-mgeneSim(genes,semData = hsGO,measure ="Rel",combine="BMA" )
        diag(semanticsim_in_cluster)<-0
        funcsimfeature<-reshape2::melt(semanticsim_in_cluster)
        colnames(funcsimfeature)<-c("TF","TG","semsim")
        funcsimfeature<-subset(funcsimfeature,TF%in%potentialTFname &TG%in%genesincluster& TF!=TG)
        matrixfuncsimfeature<-dcast(funcsimfeature,TG~TF,value.var = "semsim",fill=0,drop =F,
                                    fun.aggregate = min)
        rownames(matrixfuncsimfeature)<-matrixfuncsimfeature[,"TG"]
        matrixfuncsimfeature<-matrixfuncsimfeature[,-1]
        
        matprior[rownames(matrixfuncsimfeature),colnames(matrixfuncsimfeature)]<-as.data.frame(matrixfuncsimfeature)
        # Transforming functional similarity into weight
        listtf<-colnames(matrixfuncsimfeature)
        
        numbergenes<-length(allgenes)
        nbtf<-length(listtf)
        
        
        for (g1 in allgenes)
        {
          for (g2 in listtf)
          {
            prior<- matprior[g1,g2]
		if (prior!=0) 
		{	#print (paste0(g1,"-->",g2))
          	 matweight[g1,g2]<-ifelse(g1==g2,0,1/(prior**exponent))
		}
          }
        }
      }
    else
      if(priordata=="KD")
      {
        matprior=matrix(1,nr = numbergenes, nc = numbergenes)
        
        #exponent=1.1
        matprior<-as.data.frame(matprior)
        #test<-as.data.frame(matrixbidingprior)
        rownames(matprior)<-allgenes
        colnames(matprior)<-allgenes
        fileinfoKD<-"../data/data_human/final_data_hum_reg_network/Hela_data/KD_data/all_diff_exprs_TF.tsv"
        dataTFKD<-fread(fileinfoKD,sep="\t" ,colClasses=c("character","character","numeric","numeric","numeric")
                        ,nrows=13000000,fill=T,header=T,stringsAsFactors = F,quote = "")
        clusterdataTFKD<-subset(dataTFKD,TF%in%potentialTFname & Gene %in% genesincluster,
                                select=c("TF","Gene","adj.P.Val"))
        colnames(clusterdataTFKD)<-c("TF","TG","Pval")
        matrixKDdiffexpr<-reshape2::acast(clusterdataTFKD,TG~TF,value.var = "Pval",fill=1,drop =F,
                                    fun.aggregate = min)
        
        
        matprior[rownames(matrixKDdiffexpr),colnames(matrixKDdiffexpr)]<-as.data.frame(matrixKDdiffexpr)
        # Transforming functional similarity into weight
        listtf<-colnames(matrixKDdiffexpr)
        
        
        for (g1 in allgenes)
        {
          for (g2 in listtf)
          {
            prior=matprior[g1,g2]
            #print(pval)
            # prob=ifelse(g1==g2,0,(intergral(beta=beta,pval=pval)))
            # matweight[g1,g2]=ifelse(prob>=0.8,0.0001,1)
            matweight[g1,g2]=ifelse(g1==g2,0,1/(integral(beta=beta,pval=prior,lambdamin=lambdamin,lambdamax=lambdamax))**exponent)
          }
        }
        
        
      }
  }
print (any(is.infinite(matweight)))
  ################### Then 
  final_expr_mat<-scale(t(considered_exprdat))
  #listtf<-potentialTFname
  listgenes<-genesincluster
  matweight[matweight==0]=max(matweight[matweight!=0])
print(nbBootstrap)
if (priordata=="none")
{
print("We run BENIN without prior knowledge")
listtf<-potentialTFname
exectime<-system.time(resbenin<-applybootstrapbenin(X=final_expr_mat,nbBoobstrap=nbBootstrap,sizenetwork=sizenetwork,normalize=TRUE,nbfolds=nbfolds,
        alphaenet=alphaenet,lmean=lmean,listtf=listtf,allgenes=listgenes,parallel=TRUE))[3]
}
else
{
	print(listtf)
  exectime<-system.time(resbenin<-applybootstrapbenin(X=final_expr_mat,nbBoobstrap=nbBootstrap,matweightpk=matweight,sizenetwork=sizenetwork,normalize=TRUE,nbfolds=nbfolds,
	alphaenet=alphaenet,lmean=lmean,listtf=listtf,allgenes=listgenes,parallel=TRUE))[3]
}
print(paste0("BENIN took ",exectime," to finish the execution of the cluster",clusterid) )
#print(resbenin) 
 rescluster<-list2df(resbenin)
#filedatares=paste0("../res_in_sillico/res_human/res_hela_network/cluster_",clusterid,"/reshelanetwork_with_",priordata,"_on_cluster_",clusterid,".RData")
#save(rescluster,file=filedatares)
  rescluster$TF<-as.character(rescluster$TF)
  rescluster$TG<-as.character(rescluster$TG)
  rescluster<-subset(rescluster,TF!=TG)
  rescluster<-rescluster[order(rescluster$W,decreasing = T),]
  ## Writing res in file
  tostringalphaenet=toString(alphaenet)
  tostringalphaenet=gsub(".","",tostringalphaenet,fixed=T)
  tostringexponent=toString(exponent)
  tostringexponent=gsub(".","",tostringexponent,fixed=T)


filedatares=paste0("../res_in_sillico/res_human/res_hela_network/cluster_",clusterid,"/reshelanetwork_with_",priordata,
                 "_nbbootstrap_",nbBootstrap,"_alphaenet_",tostringalphaenet,"_exponent_",tostringexponent,"_helanetwork_",clusterid,".RData")
save(resbenin,file=filedatares)
  fileres=paste0("../res_in_sillico/res_human/res_hela_network/cluster_",clusterid,"/reshelanetwork_with_",priordata,
                 "_nbbootstrap_",nbBootstrap,"_alphaenet_",tostringalphaenet,"_exponent_",tostringexponent,"_helanetwork_",clusterid,".txt")
  print(fileres)
  write.table(rescluster,fileres,col.names = T,row.names = F,quote = F,sep = "\t")
}
