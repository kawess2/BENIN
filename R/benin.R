applycvglmnet<- function(x,y,penalty.factor=NULL,nfolds,alpha){
  if(hasArg(penalty.factor)) 
  {
   
    exectime<- system.time(fit<-cv.glmnet(x=x,y=y,
                                          family="gaussian",nfolds=nfolds,penalty.factor=penalty.factor,alpha=alpha,
                                          standardize=FALSE,intercept=FALSE,grouped=FALSE))[3]
  
  }else{
    exectime<-system.time(fit<-cv.glmnet(x=x,y=y,
                                         family="gaussian",nfolds=nfolds,alpha=alpha,
                                         standardize=FALSE,intercept=FALSE,grouped=FALSE))[3]
  }
  
  return (coef(fit, s = "lambda.min")[-1,])

  
  
}

applycvglmnetprime<- function(x,indices,y_t,penalty.factor=NULL,nfolds,alpha){
  if(hasArg(penalty.factor)) 
 {
   
     out <- tryCatch(
        { cv.glmnet(x=x,y=y_t[indices,],
                   family="gaussian",nfolds=nfolds,penalty.factor=penalty.factor,alpha=alpha,
                   standardize=FALSE,intercept=FALSE,grouped=FALSE)
}
,
        error=function(cond) {
            print("Here's the original error message after apply glmnet:")
           print(cond)
                lambda <- exp(seq(log(0.001), log(5), length.out=15))
                fit<-cv.glmnet(x=x,y=y_t[indices,],
                   family="gaussian",nfolds=nfolds,penalty.factor=penalty.factor,alpha=alpha,
                   standardize=FALSE,intercept=FALSE,grouped=FALSE, lambda=lambda)
            return(fit)
        },
	warning=function(cond) {
            #message(paste("URL caused a warning:", url))
            #print("Here's the original warning message:")
            #print(cond)
            # Choose a return value in case of warning
            #return(NULL)
        },
	finally={
        # NOTE:
})
fit<-out
  }else{
    exectime<-system.time(fit<-cv.glmnet(x=x,y=y_t[indices,],
                   family="gaussian",nfolds=nfolds,alpha=alpha,
                   standardize=FALSE,intercept=FALSE,grouped=FALSE))[3]
  }

  
  return (coef(fit, s = "lambda.min")[-1,])
  
  
}


bootstrapbeninbeta<-function (currentgene,Xscaled,nbBoobstrap=1000,matweightpk=NULL,sizenetwork,normalize=TRUE,
                              nbfolds,alphaenet,lmean=10,listtf="",parallel=TRUE)
{

  nbtimepoint<-dim(Xscaled)[1]
  reg<-unique(c(listtf,currentgene))
  considerweight<-hasArg(matweightpk)
  y_t<-as.matrix(Xscaled[(2:nbtimepoint),currentgene])
  X_t_minus_one<-as.matrix(Xscaled[(1:nbtimepoint-1),reg])
  if (considerweight)
  {

    genes_weight<-matweightpk[currentgene,reg]
     genes_weight[currentgene]<-median(genes_weight)
 genes_weight<-normalizeWeigth2(genes_weight)
  resample<-tryCatch(
	{tsboot(X_t_minus_one,R=nbBoobstrap,sim="geom",statistic = applycvglmnetprime,
                     l=lmean,y_t=y_t,penalty.factor=genes_weight,nfolds=nbfolds,alpha=alphaenet,parallel = "multicore", ncpus=2)},
	error=function(e){
                        
  print(paste0("the error here after tsboot is : ",e))
                })
 }
  else{
    resample<-tsboot(X_t_minus_one,R=nbBoobstrap,sim="geom",statistic = applycvglmnetprime,
                     l=lmean,y_t=y_t,nfolds=nbfolds,alpha=alphaenet,parallel = "multicore", ncpus=2)
  }

  
res<-tryCatch(
	{
		as.vector(colSums(resample$t!=0)/nbBoobstrap)
	},error=function(e){
			print(paste0("the error here is ",e))
		})
if (length(res)!=length(names(resample$t0)) ) 
	{
		res=vector("numeric",length(reg))
		names(res)=reg
		}
else
	{  	
	
		if (length(res)!=length(names(resample$t0)))
		{
		  warning('The length of the res vector is different than the length of names vector')
		}
		names(res)<-names(resample$t0)
		}
 return( res)
}
applybootstrapbenin<-function(X,nbBoobstrap=1000,matweightpk=NULL,sizenetwork,normalize=TRUE,nbfolds,alphaenet,lmean=10,listtf="",allgenes="",
                             parallel=TRUE,seed=225)
{
  
#   if(sizenetwork==10)
#   {
#     load("random_state_seed1001_10.RData")
#   }else
#     if(sizenetwork==100)
#     {
#     
#       load("random_state_seed1001.RData")
#     }
#   else{
# 
# load("random_state_seed1001.RData")
# 
#   }
  if(hasArg(seed))
  {
    set.seed(seed)
  }
  
  if (normalize)
  {
    Xscaled=scale(x=X,scale = apply(X, 2, sd, na.rm = TRUE))
  }
  else
  {
    Xscaled=X
  }
  considerweight<-hasArg(matweightpk)
  clus <- parallel::makeCluster(detectCores()-1,outfile="")
clusterSetRNGStream(clus,.Random.seed)
clusterEvalQ(clus, .libPaths("./library"))
  clusterExport(clus, c("detectCores","make.ends","ts.array","ts.return","isMatrix","tsboot","applycvglmnetprime","cv.glmnet","bootstrapbeninbeta","normalizeWeigth2","colSums"))  
if(considerweight)
  {
   withCallingHandlers({ output<-tryCatch({parallel::parSapply(clus,allgenes,FUN=bootstrapbeninbeta,USE.NAMES = T,simplify = F,Xscaled=Xscaled,nbBoobstrap=nbBoobstrap,
                                matweightpk=matweightpk,sizenetwork=sizenetwork,
                                nbfolds=nbfolds,alphaenet=alphaenet,lmean=lmean,listtf=listtf)},error=function(e){
                        print(paste0("the error here is ",e))
			
                })}, error=function(e) print(sys.calls()))
  }
  else
  {
    output<-parallel::parSapply(clus,allgenes,FUN=bootstrapbeninbeta,USE.NAMES = T,simplify = F,Xscaled=Xscaled,nbBoobstrap=nbBoobstrap
                                ,sizenetwork=sizenetwork,nbfolds=nbfolds,alphaenet=alphaenet,lmean=lmean,listtf=listtf)
  }
  return(output)
}
              

