#Functional analysis of the reconstructed network
options(warn=-1)
require(org.Sc.sgd.db)
jaccardIndex<-function(setReg1, setReg2)
{
  inter<-length(intersect(setReg1, setReg2))
  union<-length(union(setReg1,setReg2))
  if(union==0)
  {
    indx<-0
  }else
  {
    indx<-inter/union
  }
  
  return(indx)
}
generaterandomNetwork<-function(adjMatrix,allgenes,listtf){
  #adjMatrix<-adjMatrix[,listtf]
  row.names (adjMatrix) <- sample (row.names (adjMatrix))
  colnames(adjMatrix)<-sample(colnames(adjMatrix))
  

  for(g1 in allgenes)
  {
   #for (g2 in listtf)
    #{  
     # if (g1==g2)
      #{
        #print(g1)
        #print(g2)
        adjMatrix[g1,g1]=0
      #}
      #}
  }
  
  #print(adjMatrix[is.na(adjMatrix)])
  return(adjMatrix)
}
functionalAnalysis<-function(allgenes, global_res,listtf){
  #options(warn=-1)
  #global_res<-list2df(edgesList)
  # colnames(global_res)=c("TF","TG","W")
  global_res[(global_res[,3]>=0.5),3]<-1
  global_res[(global_res[,3]<0.5),3]<-0
  global_res[global_res[,1]==global_res[,2],3]=0
  mat_global_res<-global_res%>%dcast(TG~TF,fill=0)
  rownames(mat_global_res)<-mat_global_res[,1]
  mat_global_res<-mat_global_res[,-1]
 # mat_global_res[ is.na(mat_global_res)]<-0
  #View(mat_global_res)
    return(networkfuncAnnot(allgenes,mat_global_res,listtf))
  }
  networkfuncAnnot<-function(allgenes,adjMat,listtf)
  {
    #adjMat[is.na(adjMat)]<-0
    numbergenes<-length(allgenes)
    #numbertf<-
    matcorrgenes<-matrix(0,nr = numbergenes, nc = numbergenes)
    rownames(matcorrgenes)<-allgenes
    colnames(matcorrgenes)<-allgenes
    for (g1 in allgenes)
    {
      for (g2 in allgenes) 
      {
        if (g1!=g2)
        {
       
          listRegg1<-colnames(adjMat[g1,(adjMat[g1,]==1)])
         #print(listRegg1)
         #print("And the list of g1 is ")
       # if(g1%in%listRegg1){
       #   print("--------------------------------------------------------")
       # }
       
       # print(g2)
        listRegg2<-colnames(adjMat[g2,(adjMat[g2,]==1)])
        #print(listRegg2)
       #print(listRegg2)
       # print("And the list of g2 is ")
       # if(g2%in%listRegg2){
       #  print("-##########################################################")
       # }
       
        jaccsimm<-jaccardIndex(listRegg1,listRegg2)
        #print(jaccsimm)
        matcorrgenes[g1,g2]<-jaccsimm
       }
        
      }
    }
    matcorrgenes[matcorrgenes>0.6]=1
    matcorrgenes[matcorrgenes<=0.6]=0
    #print("######################################")
    #print(matcorrgenes)
    allcombin<-expand.grid(allgenes, allgenes,stringsAsFactors = FALSE)
    res_func_anal<-apply(allcombin,1,findAnnotationgenes,matcorrgenes )
  #View(matcorrgenes)
 
  #   matcofunc<-matrix(-1,nr = numbergenes, nc = numbergenes)
  #   rownames(matcofunc)<-allgenes
  #   colnames(matcofunc)<-allgenes
  # #View(matcorrgenes)
  #   # nbelm<-0
  #   # sumjaccindxannot<-0
  #   
  #  
  #     
  #         
  #       #print(annottableg1)
  #      # print(listannotg1)
  #         
  #      # print(annottableg2)
  #       #print(listannotg2)
  #         
  #         matcofunc[g1,g2]<-jaccindxannot
  #       
  #           nbelm<-nbelm+1
  #         sumjaccindxannot<-jaccindxannot+sumjaccindxannot
  #     
  #   
  #   }
  # #print(mean(matcofunc[matcofunc>0]))
  return(res_func_anal)
    
  }
  
  findAnnotationgenes<-function(combin, matcorrgenes)
  {
    #print(combin)
    g1<-combin[1]
    #print(gene1)
    g2<-combin[2]
    jaccindxannot<--1
    if (matcorrgenes[g1,g2]==1)
    {
      #print(g1==g2)
      annottableg1<-AnnotationDbi::select(org.Sc.sgd.db, g1, columns = "GO",ketypes="ORF")
      annottableg1<-subset(annottableg1,ONTOLOGY=="BP")
      listannotg1<-unique(annottableg1[,"GO"])
      annottableg2<-AnnotationDbi::select(org.Sc.sgd.db, g2, columns = "GO",ketypes="ORF")
      annottableg2<-subset(annottableg2,ONTOLOGY=="BP")
      listannotg2<-unique(annottableg2[,"GO"])
      jaccindxannot<-jaccardIndex(listannotg1,listannotg2)
    }
    return (jaccindxannot)
  }
  functionalAnalisisRandNetwork<-function(allgenes,listreglink=NULL,listtf,adjMatrix=NULL)
  {
    consideradjMat<-hasArg(adjMatrix)
    if (!consideradjMat)
    {
      global_res$TF<-as.character(global_res$TF)
      global_res$TG<-as.character(global_res$TG)
      global_res[(global_res[,3]>=0.5),3]<-1
      global_res[(global_res[,3]<0.5),3]<-0
      global_res[global_res[,1]==global_res[,2],3]=0
      adjMatrix<-global_res%>%dcast(TG~TF,fill=0)
      rownames( adjMatrix)<- adjMatrix[,1]
      adjMatrix<- adjMatrix[,-1]
    }
   
   #print(adjMatrix)
    foreach(i=1:1,.combine='c',.packages=c("AnnotationDbi","org.Sc.sgd.db"),.export=c("generaterandomNetwork", "networkfuncAnnot","jaccardIndex","findAnnotationgenes"))%dopar%
     #for(i in 1:100)
       {
        randnetwork<-generaterandomNetwork(adjMatrix,allgenes,listtf)
       # print(randnetwork)
       res<-networkfuncAnnot(allgenes,randnetwork,listtf)
        #names(res)<-i
        res=list(res)
        names(res)<-i
        res
        #allres<-c(allres,res)
     }
    
    #return(allres)
  }
#
   rand_net_anal_wprior<-functionalAnalisisRandNetwork(allgenes,global_res,listtf)
   lmean=lapply(rand_net_anal_wprior, function(x) mean(x[x>-1]))
   mean(unlist(lmean))
#   rand_net_anal_noprior<-functionalAnalisisRandNetwork(allgenes,global_resno_prior,listtf)
#   lmeannoprior=lapply(rand_net_anal_noprior, function(x) mean(x[x>-1]))
#   mean(unlist(lmeannoprior))
#   rand_net_anal_irafnet<-functionalAnalisisRandNetwork(allgenes,global_res_irafnet,listtf)
#   lmeanirafnet=lapply(rand_net_anal_irafnet, function(x) mean(x[x>-1]))
#   mean(unlist(lmeanirafnet))
  # funanal_irafnet<-functionalAnalysis(allgenes, global_res_irafnet,listtf)
  # mean(funanal_irafnet[funanal_irafnet>-1])
  
  # funanal_gelnet<-functionalAnalysis(allgenes, global_res_gelnet,listtf)
  # rand_net_anal_gelnet<-functionalAnalisisRandNetwork(allgenes,global_res_gelnet,listtf)
  # lmeangelnet=lapply(rand_net_anal_gelnet, function(x) mean(x[x>-1]))
  # 
  # label=lapply(1:100,FUN = function(x)(paste("rand",x)) )
  # replaceval<-function(x){
  #   #x=l[[ind]]
  #   x[x==-1]<-NA
  #   x
  # }
  # for (i in 1:10)
  # {
  #   print(testrand_net_anal_wprior[i])
  # }
  # 
  # geomean<-function(values){
  #   prod(values)**(1/length(values))
  # }
  # 
  # #testrand_net_anal_wprior<-lapply(testrand_net_anal_wprior,FUN = replaceval)
  # wilcox.test(unlist(funanal,use.names = F),unlist(rand_net_anal_wprior[[2]],use.names = F),paired=T)
  # method="benin"
  # label=lapply(1:100,FUN = function(x)(paste("rand",x)) )
  # labelmethod<-c(method,label)
  # factorlabel<-labelmethod
  # labelmethod<-rep(labelmethod,each=625)
  # labels<-factor(unlist(labelmethod,use.names = F),labels = unlist(factorlabel,use.names = F))
  # dataforwiltest<-c(unlist(funanalwprior,use.names = F),unlist(rand_net_anal_wprior,recursive = T,use.names = F))
  # testwilcox<-pairwise.wilcox.test(dataforwiltest, labels,paired = T,p.adjust="bonf")
  # exp(mean(log(testwilcox$p.value[,1])))
  # label_noprior<-c("")
  # 
  # 
  # rand_net_anal_networkBMA<-functionalAnalisisRandNetwork(allgenes,global_res_networkBMA,listtf)
  # lmeannetworkBMA=lapply(rand_net_anal_networkBMA, function(x) mean(x[x>-1]))
  # mean(unlist(lmeannetworkBMA))
  # 
  # 
