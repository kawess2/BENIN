source("bootfuns.q")
source("bootpracs.q")
source("benin2.R")
source("list2df.R")
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
library(readr)
#library(pkgconfig,lib.loc='library/')
library(GOSemSim)
library(org.Hs.eg.db)
#library(readr)

# library(Rcpp,lib.loc='library/')
# library("leaps",lib="library/")
# library("robustbase",lib="library/")
# library("inline",lib="library/")
# library(rrcov,lib="library/")
# library("BMA",lib="library/")
# library(RcppArmadillo,lib="library/")
# library(RcppEigen,lib="library/")
# library("networkBMA",lib="library/")
# library(DBI, lib = "library/")
# library(Biobase, lib = "library/")
# library(S4Vectors, lib = "library/")
# library(IRanges, lib = "library/")
# library(RSQLite, lib = "library/")
# library(AnnotationDbi, lib = "library/")
# library(org.Sc.sgd.db,lib.loc='library/')
# library('foreach',lib.loc='library/')
# library('snow',lib.loc='library/')
# library('iterators',lib.loc='library/')
# library('doSNOW',lib.loc='library/')
# library('glmnet',lib.loc='library/')
# library('parallel',lib.loc='library/')
# library('boot',lib.loc='library/')
# #library(Rcpp,lib.loc='library/')
# library(tseries,lib.loc='library/')
# library(zoo,lib.loc='library/')
# library(xts,lib.loc='library/')
# library(boot,lib.loc='library/')
# library(DMwR,lib.loc='library/')
# library(caTools,lib.loc='library/')
# library(RGeode,lib.loc='library/')
# library(crayon,lib.loc='library/')
# library(dplyr,lib.loc='library/')
# library(data.table,lib.loc='library/')
# library(reshape2,lib.loc='library/')
# library(precrec,lib.loc='library/')
# library('parallel')
library(org.Sc.sgd.db)
library(tidyverse)
library(seqinr)
#library(universalmotif)
#library(tidyr)
#library(mice)
#library(DMwR)
#library(reshape2)
#library(precrec)
#library('foreach')
#library('snow')
#library('iterators')
#library('doSNOW')
#library('glmnet')
#library('parallel')
#library('boot')
#library(doRNG)
#library(dplyr)
### Function definition
writemotif<-function(indx,allmotifsfilename,matallmotif)
{
  currentmotifID=matallmotif[indx,]$Motif_ID
  #print(currentmotifID)
  currentmotiffilename=paste0("../data/data_saccharomyces/Saccharomyces_cerevisiae_2020_04_15_10-36_pm/pwms_all_motifs/",currentmotifID,".txt")
  head=paste0("\nTF\t",paste0(matallmotif[indx,]$TF_ID,"\n"),"TF Name\t",paste0(matallmotif[indx,]$TF_Name,"\n"),"Gene\t",
              paste0("","\n"), "Motif\t",paste0(matallmotif[indx,]$Motif_ID,"\n"),"Family\t",paste0(matallmotif[indx,]$Family_Name,"\n"), 
              "Species\t",paste0(matallmotif[indx,]$TF_Species,"\n"))
  #print(currentmotiffilename)
  currentmotif=read.table(currentmotiffilename,header = F, sep = "\t",stringsAsFactors = F)
  
  #print(currentmotif)
  print(head)
  lenmotif=dim(currentmotif)[1]
  if (lenmotif>1)
  {
    write(head,allmotifsfilename,append = T)
    write.table(currentmotif,allmotifsfilename,append = T, quote=F, sep="\t", col.names = F,row.names = F)
  }
  
  
}
cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}
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

# gc()
# rm(list = ls(all.names = TRUE))
# ######### Beginning of the script
# ysexprdata<- read.csv("../data/data_saccharomyces/time_expression_data/combined.txt",sep = "\t",stringsAsFactors = F,row.names = 1,quote = "")
# listcell_cycle_gene<-read.csv("../data/data_saccharomyces/time_expression_data/cell_cycle_genes_saccharomyces.txt",sep = "\t",stringsAsFactors = F,
#                               header = T)
# symbols<-unique(listcell_cycle_gene$SGD)
# mappedgenes<-AnnotationDbi::select(org.Sc.sgd.db, as.character(symbols), c("ENSEMBL","GENENAME","SGD","ENTREZID"), "GENENAME")
# notmapped<-subset(mappedgenes, (is.na(ENSEMBL) |is.na(GENENAME)|is.na(SGD)|is.na(ENTREZID)))
# mappedgenes2<-AnnotationDbi::select(org.Sc.sgd.db, as.character(notmapped$COMMON), c("COMMON","ENSEMBL","GENENAME","SGD","ENTREZID"), "ENSEMBL")
# mappedgenes[mappedgenes$COMMON%in%mappedgenes2$ENSEMBL, "ENSEMBL"]<-unique(mappedgenes2$ENSEMBL)
# mappedgenes3<-AnnotationDbi::select(org.Sc.sgd.db, c("HDR1", "PTS1","HSN1","SDL1" , "SPH1"  , "SNR17A"), 
#                                     c("COMMON","ENSEMBL","GENENAME","SGD","ENTREZID"), "GENENAME")
# 
# mappedgenes[mappedgenes$COMMON=="HDR1","ENSEMBL"]<-"YBR138C"
# mappedgenes[mappedgenes$COMMON=="PTS1","ENSEMBL"]<-"YDR055W"
# mappedgenes[mappedgenes$COMMON=="HSN1","ENSEMBL"]<-"YHR127W"
# mappedgenes[mappedgenes$COMMON=="SDL1","ENSEMBL"]<-"YIL168W"
# mappedgenes[mappedgenes$COMMON=="SPH1","ENSEMBL"]<-"YLR313C"
# mappedgenes[mappedgenes$COMMON=="SNR17A","ENSEMBL"]<-"YOR235W"
# mappedgenes<-subset(mappedgenes,mappedgenes$COMMON!="")
# ### Saving genes mapping for latter use
# write.table(mappedgenes,"../data/data_saccharomyces/time_expression_data/cellcyclegenemapping.txt",col.names = T,row.names = F,
#             quote = F, sep="\t")
# 
# ### Subsetting the expression data to cell cycle genes 
# finalcellcyclegenes<-unique(mappedgenes$ENSEMBL[!is.na(mappedgenes$ENSEMBL)])
# 
# #filteredexprdata<-subset(subysexprdata,)
# filteredexprdata<-ysexprdata[sapply(ysexprdata, function(x) !all(is.na(x)))] 
# #filteredexprdata<-t(filteredexprdata)
# 
# # imputed_filteredexprdata <- mice(filteredexprdata, m=5, maxit = 50, method = 'pmm', seed = 500,printFlag = FALSE)
# # summary(imputed_filteredexprdata)
# # finalexprdata<-complete(imputed_filteredexprdata,2)
# cleanexprdata<-knnImputation(filteredexprdata,k=12)
# 
# 
# # We removed this genes because they are unknown in yeastract, and dubious in SGD or even deleted (YCLX09W)
# modfinalcellcyclegenes= setdiff(finalcellcyclegenes,c("YML035C-A","YCLX09W","YCL022C"))
# write.table(modfinalcellcyclegenes,"../data/data_saccharomyces/time_expression_data/finalcellcyclegenes.txt",col.names = F,row.names = F,
#             quote = F, sep="\t")
# 
# finalcleanexprdata<-cleanexprdata[modfinalcellcyclegenes,]
# 
# write.table(finalcleanexprdata,"../data/data_saccharomyces/time_expression_data/finalcleanexprdatasacc.txt",col.names = T,row.names = T,
#             quote = F, sep="\t")
# # #clustering
# exprdatfile<-"../data/data_saccharomyces/time_expression_data/finalcleanexprdatasacc.txt"
# 
# final_sacc_exprdata<-read.table(exprdatfile,header=T, stringsAsFactors =F,quote="",sep="\t",row.names =1)
# # row.names(final_annot_human_exprdata_unique_id)<-final_annot_human_exprdata_unique_id$Symbol
# # final_annot_human_exprdata_unique_id<-final_annot_human_exprdata_unique_id[,-1]
# # considered_exprdat<-final_annot_human_exprdata_unique_id[,39:86]
# # #write.table(row.names(considered_exprdat),"../data/data_human/final_data_hum_reg_network/Hela_data/considered_genes_hela_cell_cycle.txt",col.names = F,
# # #            row.names = F,quote=F)
# scaled_considered_exprdat<-scale(final_sacc_exprdata)
# corr_sacc_expr_data<- cor(t(scaled_considered_exprdat))
# gene_dist<-dist(corr_sacc_expr_data)
# gene_hclust <- hclust(gene_dist,method="average")
# plot(gene_hclust, labels = FALSE)
# nbcluster=15
# # gene_cluster <- cutree(gene_hclust, k = nbcluster) %>% 
# #   #   # turn the named vector into a tibble
# #    enframe() %>% 
# #   #   # rename some of the columns
# #     rename(gene = name, cluster = value)
# # for(idcluster in 1:nbcluster)
# # {
# #   namefilecluster<-paste0("../data/data_saccharomyces/modules/cluster",idcluster,"_avg_clust.txt")
# #   genesincluster<-subset(gene_cluster,cluster==idcluster)$gene
# #   namecluster=paste0("cluster_",idcluster)
# #   #print(length(genesincluster))
# #   #print(mean(cor(t(final_sacc_exprdata[genesincluster,]))))
# #   #write.table(genesincluster, file=namefilecluster,row.names = F,
# #   #            quote=F,sep= "\t",col.names = F)
# #   if (idcluster==1)
# #   {
# #     print("heregsgsgs")
# #     myclusterlist=list(genesincluster)
# #     names(myclusterlist)<-namecluster
# #     print(myclusterlist)
# #     allclustersdata=myclusterlist
# #   }else{
# #     myclusterlist=list(genesincluster)
# #     names(myclusterlist)<-namecluster
# #     allclustersdata= c(allclustersdata,myclusterlist)
# #   }
# # }
# 
# 
# ###### For spectral clustering
# listcluster<-specc(scaled_considered_exprdat,kernel = "rbfdot", centers= nbcluster)
# 
# for(idcluster in 1:nbcluster)
# {
#   namefilecluster<-paste0("../data/data_saccharomyces/modules/cluster",idcluster,"_specclust_clust.txt")
#   genesincluster<-names(listcluster[listcluster==idcluster])
#   namecluster=paste0("cluster_",idcluster)
#   #print(length(genesincluster))
#   #print(mean(cor(t(scaled_considered_exprdat[genesincluster,]))))
#   write.table(genesincluster, file=namefilecluster,row.names = F,
#               quote=F,sep= "\t",col.names = F)
#   if (idcluster==1)
#   {
#     #print("heregsgsgs")
#     myclusterlist=list(genesincluster)
#     names(myclusterlist)<-namecluster
#     print(myclusterlist)
#     allclustersdata=myclusterlist
#   }else{
#     myclusterlist=list(genesincluster)
#     names(myclusterlist)<-namecluster
#     allclustersdata= c(allclustersdata,myclusterlist)
#   }
# }
# #finalallclustersdata=t(plyr::ldply(allclustersdata, rbind))
# finalallclustersdata=t(plyr::ldply(allclustersdata, rbind))
# colnames(finalallclustersdata)<-finalallclustersdata[1,]
# finalallclustersdata<-finalallclustersdata[-1,]
# namefilecluster<-paste0("../data/data_saccharomyces/modules/cluster_all_spectralclust_clust.txt")
# write.table(finalallclustersdata,file=namefilecluster,row.names = F,
#            quote=F,sep= "\t",col.names = T)
# 
# # # We split the promoter sequences in different files corresponding to each clusters
# ### If a "Warning message:In readLines(file) : incomplete final line found  Navigate to the very last line of the file 
# ##### and return (press return)
#  mypromoter<-read.fasta(file = "../data/data_saccharomyces/sequences/promoter_sequences_all_Yeastract1000BP.fasta.txt",
#                         seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)
#  mappinggenefromyeastract<-read.table(file = "../data/data_saccharomyces/sequences/mappingYeastCellCycle.txt",
#                                       header = T,sep="\t",stringsAsFactors = F,quote="")
# for (id_cluster in 1:nbcluster)
#   {
#     name_file_seq_in_cluster<-paste0("../data/data_saccharomyces/sequences/cluster_",id_cluster ,".fa")
#     namefilecluster<-paste0("../data/data_saccharomyces/modules/cluster",id_cluster,"_specclust_clust.txt")
#     #print(namefilecluster)
#     genesincluster<-read.table(namefilecluster, quote = "", sep="\t",
#                                header = F,stringsAsFactors = F)
#     #print((genesincluster$V1))
#     namesseq<-subset(mappinggenefromyeastract,ORF.Name%in%genesincluster$V1)$Gene.Name
#     sequences_in_files<-mypromoter[namesseq]
#     #print(sequences_in_files)
#     seqname<-names(sequences_in_files)
#     #print(seqname)
#     submappinmappinggenefromyeastract<-subset(mappinggenefromyeastract,Gene.Name%in%namesseq) 
#     #submappinmappinggenefromyeastract<-submappinmappinggenefromyeastract[
#      # with(submappinmappinggenefromyeastract, order(namesseq)),
#     #  ]
#     #View(submappinmappinggenefromyeastract)
#     namesseq<-submappinmappinggenefromyeastract$ORF.Name
#     #print(length(namesseq))
#     write.fasta(sequences_in_files,namesseq,name_file_seq_in_cluster)
# }
#  #We want to know the potential TF:
#  KOTargetGenes=read.table("../data/data_saccharomyces/KOdata/gkq232_SuppTable1B_KOTargetGenes_Matrix_PValue.dat.txt", sep=" ")
#  symbols<-toupper(rownames( KOTargetGenes))
#  mappedTF<-AnnotationDbi::select(org.Sc.sgd.db, as.character( symbols), c("ENSEMBL","GENENAME","SGD","ENTREZID"), "GENENAME")
#  mappedTF[mappedTF$GENENAME=="RCS1","ENSEMBL"]="YGL071W"
# mappedTF[mappedTF$GENENAME=="CAF17","ENSEMBL"]="YJR122W"
# mappedTF[mappedTF$GENENAME=="RLR1","ENSEMBL"]="YNL139C"
#  mappedTF[mappedTF$GENENAME=="FLO8","ENSEMBL"]="YER109C"
# mappedTF[mappedTF$GENENAME=="RIS1","ENSEMBL"]="YOR191W" 
# mappedTF[mappedTF$GENENAME=="ZMS1","ENSEMBL"]="YJR127C"
# mappedTF[is.na(mappedTF$ENSEMBL),"ENSEMBL"]<-mappedTF[is.na(mappedTF$ENSEMBL),"GENENAME"]
# ### Use Yeasttract to convert ID it is better and complete
# locationdata<-read.csv("../data/data_saccharomyces/location_data_saccarhomyces/binding_by_gene.tsv",sep = "\t",stringsAsFactors = F)
# locationdata2<-read.csv("../data/data_saccharomyces/location_data_saccarhomyces/other_location_data_sacc_cere/pvalbygene_forpaper_abbr.txt"
#                         ,sep = "\t",stringsAsFactors = F)
# write.table(colnames(locationdata2)[4:207],"../data/data_saccharomyces/time_expression_data/TFsfromlocdata2.txt",col.names = F,row.names = F,
#             quote = F, sep="\t")
# write.table(unlist(as.list(unname(locationdata[1,!is.na(locationdata[1,])])))[5:117],
#             "../data/data_saccharomyces/time_expression_data/TFsfromlocdata1.txt",col.names = F,row.names = F,
#             quote = F, sep="\t")
# write.table(toupper(rownames( KOTargetGenes)),
#             "../data/data_saccharomyces/time_expression_data/TFsfromKOdata.txt",col.names = F,row.names = F,
#             quote = F, sep="\t")
# 
# 
# # I  copy paste the name from the second locatioon data into the file of unmapped genes so be carefull
# #mappedTFfromlocdata<-AnnotationDbi::select(org.Sc.sgd.db, as.character( locationdata$orf), c("ENSEMBL","GENENAME","SGD","ENTREZID"), "GENENAME")
# #symbols<-subset(mappedTFfromlocdata,is.na(ENSEMBL),select=c("GENENAME"))$GENENAME
# #mapp<-AnnotationDbi::select(org.Sc.sgd.db, as.character( symbols), c("ENSEMBL","GENENAME","SGD","ENTREZID"), "ENSEMBL")
# 
# allnotmappedgenes<-union(union(union(toupper(rownames( KOTargetGenes)),toupper(locationdata$orf)),
#                          unlist(as.list(unname(locationdata[1,!is.na(locationdata[1,])])))),colnames(locationdata2)[4:207])
# write.table(allnotmappedgenes,"../data/data_saccharomyces/time_expression_data/allgenesnotmapped.txt",col.names = F,row.names = F,
#             quote = F, sep="\t")
# #mappingallTF<-read.table("../data/data_saccharomyces/time_expression_data/mapping_all_TF.txt",header = T,
# #      quote = "", sep="\t",stringsAsFactors = F)mappsomegenesfromyeasttract
# TFinlocdata1<-unlist(as.list(unname(locationdata[1,!is.na(locationdata[1,])])))
# TFinlocdata2<-colnames(locationdata2)[4:207]
# mappingallTF<-read.table("../data/data_saccharomyces/time_expression_data/mappsomegenesfromyeasttract.txt",header = T,
#      quote = "", sep="\t",stringsAsFactors = F)
# #mappedorf=unique(subset(mappingallTF, mappingallTF$InsertedList%in%locationdata$orf,select="ORFName")$ORFName)
# #mappedorf1<-subset(mappingallTF,mappingallTF$Inserted.List%in%TFinlocdata1,select="ORF.Name")$ORF.Name
# mappedorf2<-subset(mappingallTF,mappingallTF$Inserted.List%in%TFinlocdata2,select="ORF.Name")$ORF.Name
# #potentialTFname1<-intersect(mappedorf1,mappedTF$ENSEMBL)
# potentialTFname2<-intersect(mappedorf2,mappedTF$ENSEMBL)
# #potentialTFname<-intersect(union(mappedorf2,mappedorf1),mappedTF$ENSEMBL)
# mappingTFcellcycle<-subset(mappingallTF,ORFName%in%potentialTFname)
# 
# write.table(mappingTFcellcycle,"../data/data_saccharomyces/time_expression_data/mappingTFcellcycle.txt",col.names = T,row.names = F,
#             quote = F, sep="\t")
# filepotentialTF<-paste0("../data/data_saccharomyces/time_expression_data/sacc_TF_list.txt")
# write.table(potentialTFname,filepotentialTF,col.names = T,row.names = F,
#             quote = F, sep="\t")
# 




# # ## Getting the list of transcription factors
# fileTFInfomotif<-paste0("../data/data_saccharomyces/Saccharomyces_cerevisiae_2020_04_15_10-36_pm/TF_Information.txt")
# TFInfomotif<-read.table(fileTFInfomotif, header = T, stringsAsFactors = F, sep="\t",fill = T)
# mappingTFCISBP<-read.table("../data/data_saccharomyces/Saccharomyces_cerevisiae_2020_04_15_10-36_pm/cisBPTFidMapping.csv",
#                            header = T, stringsAsFactors = F, sep=" ",fill = T)
# allinfoTFInfomotif<-merge(TFInfomotif,mappingTFCISBP,by.x="TF_Name",by.y="InsertedList",all.x=T)
# 
# suballinfoTFInfomotif<-subset(allinfoTFInfomotif,allinfoTFInfomotif$ORFName%in%modfinalcellcyclegenes)
# #write.table(suballinfoTFInfomotif$ORFName,"../data/data_saccharomyces/Yeast_cell_cycle_TF_havingmotif.txt",col.names = F, row.names = F,quote=F)
# 
# # 2) First we get the list of motifs file
# mymotifsfilelist=list.files("../data/data_saccharomyces/Saccharomyces_cerevisiae_2020_04_15_10-36_pm/pwms_all_motifs",no.. = T)
# listmotifid=lapply(mymotifsfilelist,FUN=function(x){gsub('.txt','',x)})
# #3 ) Get the the description of of the list of motif in the directory and which are from Transfac, chip-seq, chip-chip
# #sub_description_motif=subset(description_motif,(Motif_ID%in%listmotifid & Motif_Type %in%c("Misc","Transfac","ChIP-seq","ChIP-chip")))
# sub_description_motif<-subset(suballinfoTFInfomotif,Motif_ID !=".")
# #sub_description_motif<-subset(sub_description_motif,sub_description_motif$TF_Name %in%potentialTFname)
# #We finaly get #numbers of TF
# nbtfwmotif<-length(unique (sub_description_motif$TF_Name))
# # For a total of #motif
# nbmotif<-dim(sub_description_motif)[1]
# #We need to transform the motif in meme format
# allmotifsfilename="../data/data_saccharomyces/Saccharomyces_cerevisiae_2020_04_15_10-36_pm/allmofifscisbp.txt"
# lapply(seq(1,nbmotif),writemotif, allmotifsfilename=allmotifsfilename,matallmotif=sub_description_motif)
# cisbpmotifs<-read_cisbp(allmotifsfilename)
# memecisbpmotifsfilename="../data/data_saccharomyces/Saccharomyces_cerevisiae_2020_04_15_10-36_pm/Sacc_cere.meme"
# write_meme(cisbpmotifs,memecisbpmotifsfilename)
# # 
networkinferencepercluster<-function(priordata,clusterid,lambdamin=1,lambdamax=1000,beta=0.5,exponent=1.1,alphaenet=0.6,nbBootstrap=1000,ont="BP",methodfuncsim="Rel",sizenetwork=0)
{
  
  # We first get the data from the cluster
  clusterid= strsplit(clusterid,"_")[[1]][2]
  file_cell_cycle_genes<-"../data/data_saccharomyces/time_expression_data/finalcellcyclegenes.txt"
  cell_cycles_genes<-(read.table(file_cell_cycle_genes,header = F,quote="",sep="\t",stringsAsFactors = F))$V1
  mappingnamefile<-"../data/data_saccharomyces/sequences/mappingYeastCellCycle.txt"
  mapping_entrez_name2<-read.table(mappingnamefile,header = T,quote="",sep="\t",stringsAsFactors = F,fill=T,
                                   colClasses = c("character","character","character","character"))
  filecluster=paste0("../data/data_saccharomyces/modules/cluster_",clusterid,"_specclust_clust.txt")
  #filecluster=paste0("../data/data_human/final_data_hum_reg_network/Hela_data/modules/",clusterid,".txt")
  genesincluster<-as.vector(read.table(filecluster,stringsAsFactors = F,header = F,quote="")$V1)
  #allgenes<-unique(c(potentialTFname,genesincluster))
  
  mappingallpotentialgenes<-read.table("../data/data_saccharomyces/mappingallpotentialgenesused.txt",sep="\t",
                                       header = T, stringsAsFactors = F, quote = "")
  #Reading expression data 
  exprdatfile<-"../data/data_saccharomyces/time_expression_data/finalcleanexprdatasacc.txt"
  final_annot_yeast_exprdata_unique_id<-read.table(exprdatfile,header=T, stringsAsFactors =F,
                                                   fill=T,quote="",sep="\t",row.names = 1)
  #row.names(final_annot_human_exprdata_unique_id)<-final_annot_human_exprdata_unique_id$Symbol
  #final_annot_human_exprdata_unique_id<-final_annot_human_exprdata_unique_id[,-1]
  considered_exprdat<-final_annot_yeast_exprdata_unique_id
  #filenamegs<-"../data/data_human/final_data_hum_reg_network/Hela_data/final_human_goldstandard.txt"
  #human_gs<-read.table(filenamegs,header = T,stringsAsFactors = F,sep="\t")
  #Now we deal with each prior 
  #matpval<-matrixbidingprior
  
  matweight=matrix()
  
  if (priordata=="TFBS")
  {
    file_potential_tf<-"../data/data_saccharomyces/Yeast_cell_cycle_TF_havingmotif.txt"
    potentialTFname<-as.vector(read.table(file_potential_tf,stringsAsFactors = F,quote = "",header = F)$V1)
    allgenes<-union(potentialTFname,genesincluster)
    #analysing FIMO output
    numbergenes<-length(allgenes)
    matweight=matrix(0,nr = numbergenes, nc = numbergenes)
    rownames(matweight)<-allgenes
    colnames(matweight)<-allgenes
    matprior=matrix(1,nr = numbergenes, nc = numbergenes)
    
    #exponent=1.1
    matprior<-as.data.frame(matprior)
    #test<-as.data.frame(matrixbidingprior)
    rownames(matprior)<-allgenes
    colnames(matprior)<-allgenes
    resfilename=paste0("grep -v '^#' ../data/data_saccharomyces/sequences/res_promoter_scanning/cluster_",clusterid,"/fimo.tsv")
    resfimo<-fread(resfilename, header = T, stringsAsFactors = F,sep = "\t",colClasses = c("character","character","character","character","numeric",
                                                                                           "character","numeric","numeric","numeric","character"),
                   blank.lines.skip=T)
    fileTFInfomotif<-paste0("../data/data_saccharomyces/Saccharomyces_cerevisiae_2020_04_15_10-36_pm/",
                            "allTF_Information.txt")
    TFInfomotif<-read.table(fileTFInfomotif, header = T, stringsAsFactors = F, sep="\t",fill = T,comment.char = "")
    sub_description_motif<-subset(TFInfomotif,Motif_ID !=".")
    sub_description_motif<-subset(sub_description_motif,sub_description_motif$ORFName %in%potentialTFname)
    resfimo<-as.data.table(resfimo)
    sub_description_motif<-as.data.table(sub_description_motif)
    #resfimo$TGENSEMBL<-unlist(lapply(resfimo$sequence_name, FUN=function(x)strsplit(x,"_")[[1]][3]))
    resfimo$TGENSEMBL<-resfimo$sequence_name
    #setkey(resfimo, motif_id)
    #setkey(sub_description_motif, TF_ID)
    #df3 <- resfimo[sub_description_motif, nomatch = 0,allow.cartesian=TRUE]
    #print("fhfjfjfjffkfkfklfklfkf")    
    final_res_fimo<-inner_join(resfimo,sub_description_motif[,c(2,29,30)],by=c("motif_id"="TF_ID"))
    unique(final_res_fimo$TGENSEMBL)
    unique (final_res_fimo$ORFName)
    #mapping_genes_from_uscs<-read.table("../data/data_human/final_data_hum_reg_network/Hela_data/sequences/mapping_gene_uscs_genome.txt",header = F,
    #                                    stringsAsFactors = F,fill = T,sep="\t")
    #sub_mappingallpotentialgenes<-unique(subset(mappingallpotentialgenes,mappingallpotentialgenes$Inserted.List%in%cell_cycles_genes
     #             ,select= c("Inserted.List", "ORF.Name")))
    
    
    #colnames(sub_mappingallpotentialgenes)<-c("TGENSEMBL","REFSEQTG","TG_Name","REFSEQTG2")
    #gc(resfimo)
    #gc(mapping_genes_from_uscs)
    #final_res_fimo<-inner_join(final_res_fimo,sub_mapping_genes_from_uscs[,c("TGENSEMBL","TG_Name")],by="TGENSEMBL")
    #print(colnames(bindingfeature))  
    bindingfeature<-final_res_fimo[,c("ORFName","TGENSEMBL","q-value")]
    colnames(bindingfeature)<-c("TF_Name","TG_Name","qvalue")
    print(colnames(bindingfeature)) 
    #colnames(bindingfeature)<-c("TF_Name","DBID","TG_Name","TGENSEMBL","DBID","qvalue")
    #bindingfeature<-bindingfeature[order(bindingfeature$qvalue,decreasing = F),]
    #bindingfeature<-bindingfeature[,c("TF_Name","TG_Name","qvalue")]
    #print("fhfjfjfjffkfkfklfklfkf")
    
    matrixbidingprior<-reshape2::dcast(bindingfeature,TG_Name~TF_Name,value.var = "qvalue",fill=1,drop =F,fun.aggregate = min)
    rownames(matrixbidingprior)<-matrixbidingprior[,"TG_Name"]
    matrixbidingprior<-matrixbidingprior[,-1]
    #matpval<-matrixbidingprior
    
    #test<-as.data.frame(matrixbidingprior)
    
    matprior[rownames(matrixbidingprior),colnames(matrixbidingprior)]<-as.data.frame(matrixbidingprior)
    # Transforming location data into probabilities
    listtf<-colnames(matrixbidingprior)
    
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
    
  }else{
    if (priordata=="location")
    {
      locationdata2<-read.csv("../data/data_saccharomyces/location_data_saccarhomyces/other_location_data_sacc_cere/pvalbygene_forpaper_abbr.txt"
                              ,sep = "\t",stringsAsFactors = F)
      
      locationdata2<-locationdata2[,c(1,4:207)]
      locationdata2<-subset(locationdata2,ORFTG!="")
      rownames(locationdata2)<-locationdata2$ORFTG
      locationdata2<-locationdata2[,-1]
      TFinlocdata<-unique(subset(mappingallpotentialgenes,mappingallpotentialgenes$Inserted.List%in%colnames(locationdata2)
                         ,select= c("Inserted.List", "ORF.Name")))
      #TFinlocdata$ORF.Name.x<-TFinlocdata$ORF.Name
      potentialTFname<-intersect(TFinlocdata$ORF.Name,cell_cycles_genes)
      allgenes<-union(potentialTFname,genesincluster)
      #analysing FIMO output
      numbergenes<-length(allgenes)
      matweight=matrix(0,nr = numbergenes, nc = numbergenes)
      rownames(matweight)<-allgenes
      colnames(matweight)<-allgenes
      matprior=matrix(nr = numbergenes, nc = numbergenes)
      
      #exponent=1.1
      matprior<-as.data.frame(matprior)
      #test<-as.data.frame(matrixbidingprior)
      rownames(matprior)<-allgenes
      colnames(matprior)<-allgenes
      #matprior=matrix(1,nr = numbergenes, nc = numbergenes)
      
      row.names(TFinlocdata)<-TFinlocdata$Inserted.List
      finallocationdata<-t(locationdata2)
      #finallocationdata<-finallocationdata
      finallocationdata<-merge(finallocationdata,TFinlocdata,by=0,all=T,drop=F)
      row.names(finallocationdata)<-finallocationdata$ORF.Name
      
      drops <- c("ORF.Name","Row.names")
      finallocationdata<-finallocationdata[ , !(names(finallocationdata) %in% drops)]
      #finallocationdata<-finallocationdata[,-1]
      
      finallocationdata<-t(as.matrix(finallocationdata))
      finallocationdata<-subset(finallocationdata,rownames(finallocationdata)%in%allgenes)
      finallocationdata<-finallocationdata[,potentialTFname]
      storage.mode(finallocationdata)<-"numeric"
      matprior[rownames(finallocationdata),colnames(finallocationdata)]<-finallocationdata
      storage.mode(finallocationdata)<-"numeric"
      matprior[rownames(finallocationdata),colnames(finallocationdata)]<-as.data.frame(finallocationdata)
      #storage.mode(matprior)<-"numeric"
      #listtf<-colnames(matrixKDdiffexpr)
      listtf<-potentialTFname
      
      for (g1 in allgenes)
      {
        for (g2 in listtf)
        {
         #print(paste0(g1,"---",g2))
          prior=matprior[g1,g2]
          #print(class(prior))
          # prob=ifelse(g1==g2,0,(intergral(beta=beta,pval=pval)))
          # matweight[g1,g2]=ifelse(prob>=0.8,0.0001,1)
          if (!is.na(prior))
          {
            #print(prior)
            matweight[g1,g2]=ifelse(g1==g2,0,1/(integral(beta=beta,pval=prior,lambdamin=lambdamin,lambdamax=lambdamax))**exponent)
          }else
          {
            matweight[g1,g2]=0
          }
         
        }
      }
      
      
    }else
    {
      if (priordata=="KO")
      {
        KOdata =read.table("../data/data_saccharomyces/KOdata/gkq232_SuppTable1B_KOTargetGenes_Matrix_PValue.dat.txt", sep=" ")
        row.names(KOdata)<-toupper(rownames( KOdata))
        TFinKOdata<-unique(subset(mappingallpotentialgenes,mappingallpotentialgenes$Inserted.List%in%row.names( KOdata)
                        ,select= c("Inserted.List", "ORF.Name")))
        potentialTFname<-intersect(TFinKOdata$ORF.Name,cell_cycles_genes)
        allgenes<-union(potentialTFname,genesincluster)
        #analysing FIMO output
        numbergenes<-length(allgenes)
        matweight=matrix(0,nr = numbergenes, nc = numbergenes)
        rownames(matweight)<-allgenes
        colnames(matweight)<-allgenes
        matprior=matrix(nr = numbergenes, nc = numbergenes)
        
        #exponent=1.1
        #matprior<-as.data.frame(matprior)
        #test<-as.data.frame(matrixbidingprior)
        rownames(matprior)<-allgenes
        colnames(matprior)<-allgenes
        row.names(TFinKOdata)<-TFinKOdata$Inserted.List
        #finallocationdata<-t(locationdata2)
        #finallocationdata<-finallocationdata
        finalKOdata<-merge(KOdata ,TFinKOdata,by=0,all=T,drop=F)
        
        row.names(finalKOdata)<-finalKOdata$ORF.Name
        
        drops <- c("ORF.Name","Row.names")
        finalKOdata<-finalKOdata[ , !(names(finalKOdata) %in% drops)]
        #finallocationdata<-finallocationdata[,-1]
        finalKOdata<-t(as.matrix(finalKOdata))
        finalKOdata<-subset(finalKOdata,rownames(finalKOdata)%in%allgenes)
        finalKOdata<-finalKOdata[,potentialTFname]
        storage.mode(finalKOdata)<-"numeric"
        matprior[rownames(finalKOdata),colnames(finalKOdata)]<-finalKOdata
        #matprior<-
        
        
        listtf<-potentialTFname
        for (g1 in allgenes)
        {
          for (g2 in listtf)
          {
            #print(paste0(g1,"---",g2))
            prior=matprior[g1,g2]
           # print(class(prior))
            # prob=ifelse(g1==g2,0,(intergral(beta=beta,pval=pval)))
            # matweight[g1,g2]=ifelse(prob>=0.8,0.0001,1)
            if (!is.na(prior))
            {
              #print(prior)
              matweight[g1,g2]=ifelse(g1==g2,0,1/(integral(beta=beta,pval=prior,lambdamin=lambdamin,lambdamax=lambdamax))**exponent)
            }else
            {
              matweight[g1,g2]=0
            }
            
          }
        }
        
      }
    }
    if (priordata=="functional")
    {
      namefiletf<-"../data/data_saccharomyces/yeast_tf_names.tsv"
      Tfinsacc<-read.table(namefiletf,header = F,stringsAsFactors = F, sep="\t",quote="")
      potentialTFname<-intersect(Tfinsacc$V1,cell_cycles_genes)
      allgenes<-union(potentialTFname,genesincluster)
      #analysing FIMO output
      numbergenes<-length(allgenes)
      matweight=matrix(0,nr = numbergenes, nc = numbergenes)
      rownames(matweight)<-allgenes
      colnames(matweight)<-allgenes
      matprior=matrix(nr = numbergenes, nc = numbergenes)
      
      #exponent=1.1
      #matprior<-as.data.frame(matprior)
      #test<-as.data.frame(matrixbidingprior)
      rownames(matprior)<-allgenes
      colnames(matprior)<-allgenes
      matprior=matrix(0,nr = numbergenes, nc = numbergenes)
      
      #exponent=1.1
      #matprior<-matprior
      #test<-as.data.frame(matrixbidingprior)
      rownames(matprior)<-allgenes
      colnames(matprior)<-allgenes
      hsGO<-godata("org.Sc.sgd.db",keytype="ENSEMBL",ont = ont,computeIC = TRUE)
      genes<-unique(allgenes)
     print (paste0("the measure we are using is ",methodfuncsim))
	 semanticsim_in_cluster<-mgeneSim(genes,semData = hsGO,measure =methodfuncsim,combine="BMA" )
      diag(semanticsim_in_cluster)<-0
      funcsimfeature<-reshape2::melt(semanticsim_in_cluster)
      colnames(funcsimfeature)<-c("TF","TG","semsim")
      funcsimfeature<-subset(funcsimfeature,TF%in%potentialTFname &TG%in%genesincluster& TF!=TG)
      matrixfuncsimfeature<-dcast(funcsimfeature,TG~TF,value.var = "semsim",fill=0,drop =F,
                                  fun.aggregate = min)
      rownames(matrixfuncsimfeature)<-matrixfuncsimfeature[,"TG"]
      matrixfuncsimfeature<-matrixfuncsimfeature[,-1]
      
      matprior[rownames(matrixfuncsimfeature),colnames(matrixfuncsimfeature)]<-as.matrix(matrixfuncsimfeature)
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
          {	
           # print (paste0(g1,"-->",g2))
            matweight[g1,g2]<-ifelse(g1==g2,0,1/(prior**exponent))
          }
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
  print (any(matweight==0))
  print(nbBootstrap)
  
  resbenin<-applybootstrapbenin(X=final_expr_mat,nbBoobstrap=nbBootstrap,matweightpk=matweight,sizenetwork=sizenetwork,normalize=TRUE,nbfolds=10,alphaenet=0.6,lmean=10,listtf=listtf,allgenes=listgenes,
                                parallel=TRUE)
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
  
  
  filedatares=paste0("../res_in_sillico/res_saccharomyce/res_sacc_network/cluster_",clusterid,"/ressaccnetwork_with_",priordata,
                     "_nbbootstrap_",nbBootstrap,"_alphaenet_",tostringalphaenet,"_exponent_",tostringexponent,"_saccnetwork_",clusterid,".RData")
  save(resbenin,file=filedatares)
  fileres=paste0("../res_in_sillico/res_saccharomyce/res_sacc_network/cluster_",clusterid,"/ressaccnetwork_with_",priordata,
                 "_nbbootstrap_",nbBootstrap,"_alphaenet_",tostringalphaenet,"_exponent_",tostringexponent,"_saccnetwork_",clusterid,".txt")
  print(fileres)
  write.table(rescluster,fileres,col.names = T,row.names = F,quote = F,sep = "\t")
  
}










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
#### First we want to filter the expression data to the cell cycle genes

#Sgdgenename <- org.Sc.sgdGENENAME
# Get the gene names that are mapped to an ORF identifier
#mapped_genes <- mappedkeys(Sgdgenename)
# Convert to a list
#list_mapped_genes <- as.list(Sgdgenename [mapped_genes])
##Reading location data file to filter genes
# locationdata<-read.csv("../data/location_data_saccarhomyces/binding_by_gene.tsv",sep = "\t",stringsAsFactors = F)
# ### We match the set of gene that have previously been identified as cell cycle genes
# listcell_cycle_gene<-read.csv("../data/location_data_saccarhomyces/list_genes_cell_cycle.txt",sep = "\t",stringsAsFactors = F,header = F)
# #rawtgs<-c("ALG7","CDC20","CDC21","CDC5","CDC6","CLB2","CLB5","CLN1","CLN2","CTS1","EGT2","FAR1","HTA1","PCL2","SIC1")
# rawtgs<-intersect(listcell_cycle_gene$V1,locationdata$orf)
# 
# rawtfs<-c("ACE2","ASH1","FKH1","MBP1","MCM1","NDD1","STB1","SWI4","SWI5","SWI6")
# rawallgenes<-c(rawtgs,rawtfs)
# #match_tg<-list_mapped_genes[which(list_mapped_genes%in%c("ALG7","CDC20","CDC21","CDC5","CDC6","CLB2","CLB5","CLN1","CLN2","CTS1","EGT2","FAR1","HTA1","PCL2","SIC1"))]
# match_tg<-list_mapped_genes[which(list_mapped_genes%in%rawtgs)]
# listtg<-names(match_tg)
# rawtgs<-rawtgs[rawtgs%in%list_mapped_genes]
# ##### We match the set of transcription factor that have previously been identified as cell cycle genes
# match_tf<-list_mapped_genes[which(list_mapped_genes%in%c("ACE2","ASH1","FKH1","MBP1","MCM1","NDD1","STB1","SWI4","SWI5","SWI6"))]
# #listtf<-names(match_tf)
# 
# allgenes_matching<-c(match_tg,match_tf)
# ##### we download the expression data matrix
# #listtfnames<-c("ACE2","ASH1","FKH1","MBP1","MCM1","NDD1","STB1","SWI4","SWI5","SWI6")
# listtf<-unlist(lapply(listtfnames,FUN = function(x,db)names(db[(db==x)]),list_mapped_genes))
# allgenes<-union(listtg,listtf)
# ysexprdata<- read.csv("../data/time_expression_data/combined.txt",sep = "\t",stringsAsFactors = F)
# subysexprdata<-ysexprdata[(ysexprdata$X %in% allgenes),]
# filteredexprdata<-subysexprdata[sapply(subysexprdata, function(x) !all(is.na(x)))] 
# row.names(filteredexprdata )<-filteredexprdata [,1]
# filteredexprdata<-filteredexprdata [,-1]
# filteredexprdata<-t(filteredexprdata)
# # imputed_filteredexprdata <- mice(filteredexprdata, m=5, maxit = 50, method = 'pmm', seed = 500,printFlag = FALSE)
# # summary(imputed_filteredexprdata)
# # finalexprdata<-complete(imputed_filteredexprdata,2)
# cleanexprdata<-knnImputation(filteredexprdata,k=5)
# cleanexprdata<-scale(x=cleanexprdata,scale = apply(cleanexprdata, 2, sd, na.rm = TRUE))
# allgenes<-colnames(cleanexprdata)
# networkinferenceperclustersacc<-function(priordata,clusterid,lambdamin=1,lambdamax=1000,beta=0.5,exponent=1.1,alphaenet=0.6,nbBootstrap=1000,ont="BP",methodfuncsim="Rel",sizenetwork=0)
# {
#   
# }
