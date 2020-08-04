\name{launcher_fMM}
\title{A user frendly interface to fMMalgo}
\alias{launcher_fMM}
\usage{
launcher_fMM(G, data, Wweight_SR=NULL, distband.SR=NA, Wdist_LR=NULL, distband.LR=NA,
    var_reg, family=gaussian(link="identity"), prior_theta="kmeans", posterior_proba=NULL,
    var_group=NULL, prior_prevalence=FALSE,
    test.GR=FALSE, seed=NULL,
    test.ICM=FALSE, rho_ICM="auto", G_ICM=1:G, update_rho=FALSE,
    test.ICMregional=FALSE, coords=NULL, distance_ref=NULL, multiV=FALSE,
    epsilon=10^(-3), epsilon_corrSpat=epsilon*10, iter_max=100, trace_iter=FALSE)		
}
\description{
A more user frendly interface to the segmentation algorithm despite more limited options compared to \code{\link{fMMalgo}}
}
\arguments{
  \item{G}{The number of groups to be considered for the mixture. \code{integer}.}
  \item{data}{A dataframe containing the data. \code{data.frame[n,p]}.}
  \item{Wweight_SR}{The local neighborhood matrix. \code{dgCMatrix[n,n]}. Should be normalized by row (i.e. \code{rowSums(Wweight_SR)=1}). Only used if \code{test.GR} or \code{test.ICM} is \code{true}.}
  \item{distband.SR}{distband for automatic computation of \code{Wdist_SR}. \code{numeric}. Ignored if \code{Wweight_SR} is specified.}
  \item{Wdist_LR}{The regional neighborhood matrix. \code{dgCMatrix[n,n]}. Should contain the distances between the observations (\code{0} indicating infinite distance). Only used if both \code{test.ICM} and \code{test.ICMregional} are \code{true}.}
  \item{distband.LR}{distband for automatic computation of \code{Wdist_LR}. \code{numeric}. Ignored if \code{Wdist_LR} is specified.}
  
  \item{var_reg}{A vector indicating the names of the response variables. \code{character vector[g]}.}
  \item{family}{A list of \code{G} families indicating the distribution of each response variable. \code{list[[family]]}.}
  \item{prior_theta}{Initialisation of the means of each group means by random sampling (\code{NULL}), kmeans (\code{"kmeans"}) or user defined values (\code{list[[M]][G]}).}
  \item{posterior_proba}{Initialisation of the posterior membership probabilities to user defined values (\code{matrix[n,g]}). Ignored if \code{NULL}.}
 
  \item{var_group}{A vector indicating the names of the concomittant variables. \code{character vector[g]}.}
  \item{prior_prevalence}{Should a prior based on the prevalence of each group be used ? \code{logical}.}
 
  \item{test.GR}{Should Growing Region algorithm be performed ? \code{logical}.}
  \item{seed}{An indicator of the seeds to be used in the GR algorithm. \code{logical vector[n]}. Can alternatively contains the name of a column in data (\code{character}).}

  \item{test.ICM}{Should local/regional regularization be performed ? \code{logical}.}
  \item{rho_ICM}{Value of the spatial regularisation parameters : can be fixed a priori (\code{numeric} if test.ICMregional=F and \code{numeric vector[2]} if test.ICMregional=T) or estimated (\code{"init"})}
  \item{G_ICM}{Potential group merging for the regularization step (\code{interger vector[G]}) where each element must be inferior or equal to G.}
  \item{update_rho}{Should rho be re-estimated using the posterior probabilites at each step ? \code{logical}.}
 
  \item{test.ICMregional}{Should regional regularization be performed ? \code{logical}. \code{test.ICM} must be also T to be active.}
  \item{coords}{Coordinates (\code{matrix[n,.]}) or the name of columns in data (\code{character vector[.]}) giving the observation coordinates.}
  \item{distance_ref}{The distance defining the several neighborhood orders relatively to \code{Wdist_LR}. \code{numeric vector[.]}.}
  \item{multiV}{Should the regional potential range be specific to each spatial group ? \code{logical}. Default is \code{FALSE}.}

  \item{epsilon}{Convergence occurs when the relative variation of the log-likelihood or of the parameter values between two iterations is below \code{epsilon}. \code{numeric}.}
  \item{epsilon_corrSpat}{Spatial regularization begins when the relative variation of the log-likelihood or of the parameter values between two iterations is below \code{epsilon_corrSpat}. \code{numeric}.}
  \item{iter_max}{Maximum number of EM iterations. \code{integer}}

  \item{trace_iter}{Should estimation be displayed at each EM iteration ? \code{logical}.}
}
\value{
 A \code{list} containing :
  \item{Y}{a \code{data.frame} containing the response variables.}
  \item{Id}{a \code{data.frame} containing, if any the coordinates of the observations.}
  \item{theta}{a \code{list} containing the estimated regressors of the response models.}
  \item{sigma}{a \code{list} containing the estimated dispersion coefficient of the response models.}
  \item{reg}{a \code{list} containing the adjusted response (\code{$Y}) and concomittant (\code{$groupe}) models.}
  \item{beta}{a \code{vector} containing the estimated regressors of the concomittant models.}
  \item{rho_ICM}{if any, a \code{vector} containing the regularization parameters.}
  \item{posterior_proba}{a \code{matrix} containing the posterior probabilities of the observations.}
  \item{prior_proba}{a \code{matrix} containing the prior probabilities of the observations.}
  \item{proba_priorReg}{if any, a \code{matrix} containing the spatial prior probabilities of the observations.}
  \item{radius}{if any, a \code{vector} indicating the radius of the spatial groups of group G}
  \item{group}{a \code{list} indicating for each group the most probable index of observations that belong to the group.}
  \item{family}{a \code{list} containing the family used in the response models.}
  \item{critere_cv}{a \code{list} containing the convergence criteria.}
  \item{cv}{a \code{logical} indicating whether convergence occured.}
  \item{tps}{a \code{list} indicating the time at which the algorithm started and stopped.}  
}
\details{
This function require the following packages : \code{Matrix}, \code{nnet} and \code{spam}, \code{Rcpp}, \code{RcppArmadillo} in case of spatial regularization. 

See \code{\link{fMMalgo}} for other details.
}
\seealso{
\code{\link{fMMalgo}} for a more options, \cr
\code{\link{plotCv_fMM}} to plot the convergence criteria, \cr
\code{\link{calcW_fMM}} to compute the neighborhood matrix, \cr
\code{\link{rhoMF}} to estimate the regularization parameters on external data.
}
\examples{
data(Simulation_data,package="fMMseg")
n <- nrow(Simulation_data)

#### fMM non spatial ####
res_EMbase <- launcher_fMM(G=3,data=Simulation_data,
                            var_reg="Y_artlinear",family=gaussian(link="identity"),
                            epsilon=5*10^{-3})

#### fMM local regularization ####
 res_EMICM <-  launcher_fMM(G=3,data=Simulation_data,coords=c("i","j"),distband.SR=sqrt(2),
                            var_reg="Y_artlinear",family=gaussian(link="identity"),
                            test.ICM=TRUE,rho_ICM=c(6,0),
                            epsilon=5*10^{-3})
  
#### fMM local and regional regularization ####
system.time(
  res_EMICMlocreg <-  launcher_fMM(G=3,data=Simulation_data,coords=c("i","j"),distband.SR=sqrt(2),distband.LR=10,
                               var_reg="Y_artlinear",family=gaussian(link="identity"),
                               test.ICM=TRUE,rho_ICM=c(6,6),
                               test.ICMregional=TRUE,
                               epsilon=5*10^{-3})
)

 if(try(require(fields))){
  par(mfrow=c(2,2),mar=rep(1.5,4))
 image.plot(x=1:sqrt(n),y=1:sqrt(n),xlab="",ylab="",axes=FALSE,
            z=matrix(Simulation_data[,"Y_artlinear"],sqrt(n),sqrt(n)),
            main="Dataset")
 image.plot(x=1:sqrt(n),y=1:sqrt(n),xlab="",ylab="",axes=FALSE,
            z=matrix(res_EMbase$posterior_proba[,3],sqrt(n),sqrt(n)),
            main=expression(fMM[ns]))
 image.plot(x=1:sqrt(n),y=1:sqrt(n),xlab="",ylab="",axes=FALSE,
            z=matrix(res_EMICM$posterior_proba[,3],sqrt(n),sqrt(n)),
            main=expression(paste(fMM[Loc]," - weak",sep="")))
 image.plot(x=1:sqrt(n),y=1:sqrt(n),xlab="",ylab="",axes=FALSE,
            z=matrix(res_EMICMlocreg$posterior_proba[,3],sqrt(n),sqrt(n)),
            main=expression(fMM[LocReg]))
 }

## for other examples see the demo files : 
??fMMEvaluation # line 65 to 102 (simulated data)   
??fMMLesionSegmentation # line 55-67, 86-100 and 121-135  (real example)
}



                        
                        
                      