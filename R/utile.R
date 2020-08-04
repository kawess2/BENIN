appendG<-function(x)paste("G",x,sep="")
frequence<- function(x){length(which(x!=0))}
medianVal<-function(x){median(x)}
normalizeWeigth1<- function(x){(x-min(x))/(max(x)-min(x))}
normalizeWeigth2<- function(x){x/sum(x)} 

list2df <- function(edgeList, aks = NULL, full = FALSE) {
  onam <- names(edgeList)
  tiebreak <- NULL;
  if (!is.null(aks) && is.numeric(edgeList[!sapply(aks, is.null)][[1]])) {
    tiebreak <- as.vector(unlist(aks))
  }
  if (is.numeric(edgeList[!sapply(edgeList, is.null)][[1]])) {
    prob <- as.vector(unlist(edgeList))
    edgeList <- lapply(edgeList, names)
  }
  else prob <- NULL
  inam <- unlist(edgeList)
  edgeList <- lapply(edgeList, unlist)
  lev <- union(onam, unique(inam))
  
  if (full) {
    ll <- length(lev)
    mat <- data.frame(regulator = as.factor(rep(lev, 
                                                ll)), gene = as.factor(rep(lev, each = ll)))
    if (!is.null(prob)) {
      stop("not fixed to handle probabilities")
    }
  }
  else {  
    if (is.null(prob)) {
      mat <- data.frame(regulator = factor(inam, levels = lev), 
                        gene = factor(rep(onam, sapply(edgeList, length)), 
                                      levels = lev))
      colnames(mat) <- c( "TF", "TG" );
    }
    else {
      if ( is.null(tiebreak) ) {
        #print("check point")
        mat <- data.frame(regulator = factor(inam, levels = lev), 
                          gene = factor(rep(onam, sapply(edgeList, length)), 
                                        levels = lev), post.prob = prob);
        colnames(mat) <- c( "TF", "TG", "W" );
      }
      else {
        mat <- data.frame(regulator = factor(inam, levels = lev), 
                          gene = factor(rep(onam, sapply(edgeList, length)), 
                                        levels = lev), post.prob = prob, tiebreak = tiebreak);
        
        colnames(mat) <- c( "TF", "TG", "W", "Tiebreak" );
      }
    }
  }
  rownames(mat) <- NULL
  mat
}

writeEdges <-
  function(network, threshhold = .5, fileName= "edges.txt") {
    if (any(network[,3] > 1)) stop("probability must not be given as percentage")
    use  <- network[,3] >= threshhold
    if (!any(use)) stop("no edges have probability at or above threshhold")
    network <- network[use,,drop= FALSE]
    write.table( network, file=fileName, quote=FALSE, 
                 row.names=FALSE, col.names=FALSE, sep = "\t")
    invisible(network)
  }

randnumbergen<-function(x,lambda){
  
  ifelse(x,rexptr(1,lambda,c(0,1)),
         runif(1, min = 0, max = 1))
}
randnumbergenrow<-function(x,lambda)
{
  vapply(x,FUN =randnumbergen ,FUN.VALUE = numeric(1),lambda=lambda)
}
proba<- function(beta,pval){
  proba2<- function(lambda){
    (lambda*exp(-lambda*pval)*beta)/((lambda*exp(-lambda*pval)*beta)+
                                       ((1-exp(-lambda))*(1-beta)))
  }
  return(proba2)
}
integral<- function(beta,pval,lambdamin,lambdamax){
  return((integrate(proba(beta,pval),lower =lambdamin,upper = lambdamax)$value)/(lambdamax-lambdamin))
}
computeintegral<-function(lambda,pval,beta){
  (lambda*exp(-lambda*pval)*beta)/((lambda*exp(-lambda*pval)*beta)+((1-exp(-lambda))*(1-beta)))
}
integral2<- function(beta,pval,lambdamin,lambdamax)
{
  return((integrate(computeintegral,lower =lambdamin,upper = lambdamax,pval = pval, beta = beta )$value)/(lambdamax-lambdamin))
}
savedata<- function(data,file, colnames=TRUE, rownames=TRUE,sep="\t",quote=FALSE,append=FALSE)
{
  write.table(data,file=file,sep=sep,append = append
              ,col.names=colnames,row.names = rownames,quote=FALSE)
}

getcolor<-function(indx,lstedgcolor,lstedgcolor2,lstedgcolor3,color)
{
  currcol<-lstedgcolor[indx]
  currcol2<-lstedgcolor2[indx]
  currcol3<-lstedgcolor3[indx]
  if (is.na(currcol)&!is.na(currcol2)&!is.na(currcol3))
  {
    #blue
    return(color[26])
  }
  if (!is.na(currcol)&!is.na(currcol2)&is.na(currcol3))
  {
    #magenta
    return("black")
  }
  if (!is.na(currcol)&is.na(currcol2)&!is.na(currcol3))
  {
    #pink
    return(color[318])
  }
  if (!is.na(currcol)&!is.na(currcol2)&!is.na(currcol3))
  {
    #rouge
    return(color[334])
  }
  if (is.na(currcol)&is.na(currcol2)&!is.na(currcol3))
  {
    return(currcol3)
  }
  if (!is.na(currcol)&is.na(currcol2)&is.na(currcol3))
  {
    return(currcol)
  }
  if (is.na(currcol)&!is.na(currcol2)&is.na(currcol3))
  {
    return(currcol2)
  }
}
initPackage <- function(package, argument = NULL, tryAttach = FALSE, method){
  
  test.package <- requireNamespace(package, quietly = TRUE)
  if(test.package == FALSE){
    stop(method, " : this function ", if(!is.null(argument)){paste("with argument ", argument, " ", sep = "")}, "requires to have installed the ", package, " package to work \n")
  }
  if(tryAttach == TRUE && (paste("package:", package, sep = "") %in% search() == FALSE) ){
    try(attachNamespace(package))
  }
}

modcalcAUPRC <- function(x, y, subdivisions = 10000, performance = NULL, ci = TRUE, alpha = 0.05, 
                         method = "Kronrod", reltol = .Machine$double.eps^0.25){
  
  
  #### tests
  ## packages
  # initPackage("ROCR", method =  "calcAUPRC")	
  
  ## arguments
  if (optionsMRIaggr("checkArguments")) {
    
    validInteger(value = subdivisions, validLength = 1, min = 0, method = "calcAUPRC")
    validLogical(value = ci, validLength = 1, method = "calcAUPRC")
    if(!is.null(performance)){
      validClass(value = performance, validClass = "performance", superClasses = FALSE, method = "calcAUPRC")
      
      if( any(c("Precision", "Recall") %in% c(performance@x.name, performance@y.name) == FALSE) ){
        stop("calcAUPRC : wrong specification of \'performance\' \n", 
             "\'performance\' must contains \"Precision\" and \"Recall\" performance measures \n", 
             "measures in \'performance\' : \"", performance@x.name, "\" \"", performance@y.name, "\" \n")   
      }
    }
    validNumeric(value = alpha, validLength = 1, min = 0, max = 1, method = "calcAUPRC")
    validCharacter(value = method, validLength = 1, 
                   validValues = c("integrate", "Kronrod", "Richardson", "Clenshaw", "Simpson", "Romberg"),
                   refuse.NULL = TRUE, method = "calcAUPRC")
    if(method != "integrate"){
      #initPackage("pracma", argument =  "method != \"integrate\"", method =  "calcAUPRC")	
    }
    validNumeric(value = reltol, validLength = 1, min = 0, method = "calcAUPRC")
    
  }
  
  #### initialisation
  if(is.null(performance)){    
    performance <- ROCR::performance(ROCR::prediction(x, y), x.measure =  "rec", measure =  "prec")
  }
  
  M.PR <- cbind(performance@x.values[[1]], 
                performance@y.values[[1]])
  n <- nrow(M.PR)
  colnames(M.PR) <- c(performance@x.name,performance@y.name)
  M.PR <- M.PR[c(-1),, drop = FALSE]
  
  # deal with multiple values at the same threshold
  #M.PR <- cbind(Recall = unique(M.PR[,"Recall"]), 
  #              Precision = tapply(M.PR[,"Precision"], M.PR[,"Recall"], function(x){max(x)}))
  if(nrow(M.PR) == 1){return(M.PR[,"Precision"])} # cas avec un seuil 
  
  # conversion to a function 
  f  <- stats::approxfun(x = M.PR[,"Recall"], y = M.PR[,"Precision"])
  
  #### computation : integration
  if(method == "integrate"){
    f.int <- stats::integrate(f, min(M.PR[,"Recall"]), max(M.PR[,"Recall"]), subdivisions = subdivisions, rel.tol = reltol)$value
  }else{
    f.int <- pracma::integral(fun = f, xmin = min(M.PR[,"Recall"]), xmax = max(M.PR[,"Recall"]), 
                              method = method, reltol = reltol)
  }
  
  AUPRC <- f.int / (max(M.PR[,"Recall"]) - min(M.PR[,"Recall"]))
  
  if(ci == TRUE){
    eta <- log(AUPRC / (1 - AUPRC))
    tau <- 1 / sqrt(n * AUPRC * (1 - AUPRC))    
    
    inf_tempo <- exp(eta - stats::qnorm(p = 1 - alpha / 2) * tau)
    sup_tempo <- exp(eta + stats::qnorm(p = 1 - alpha / 2) * tau)
  }else{
    inf_tempo <- NA
    sup_tempo <- NA
  }
  
  #### export
  return(c(AUPRC = AUPRC, 
           IC_inf = inf_tempo / ( 1 + inf_tempo), 
           IC_sup = sup_tempo / ( 1 + sup_tempo)
  ))
  
}
writemotif<-function(indx,allmotifsfilename,matallmotif)
{
  currentmotifID=matallmotif[indx,]$Motif_ID
  #print(currentmotifID)
  currentmotiffilename=paste0("../data/data_human/Homo_sapiens_2019_05_23_10-17_am/pwms_all_motifs/",currentmotifID,".txt")
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

