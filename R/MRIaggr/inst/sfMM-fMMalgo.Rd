\name{fMMalgo}
\title{finite Mixture Model with spatial regularization}
\alias{fMMalgo}
\usage{
fMMalgo(M, G, data, Id=NULL, Wweight_SR=NULL, Wdist_LR=NULL, posterior_proba=NULL,
         formula_reg, offset_reg=NULL, family=gaussian(link="identity"), prior_theta="kmeans", prior_sigma=NULL,
         formula_group=NULL, prior_prevalence=FALSE, prior_proba=NULL,
         test.GR=FALSE, sigma_GR="auto", proba_GR=NULL, proba_GRseed=NULL, seed=NULL, n.order=3,
         test.ICM=FALSE, rho_ICM="init", G_ICM=1:G, prior_prevalenceICM=TRUE, rho_max=10, update_rho=FALSE,                            
         test.ICMregional=FALSE, coords=NULL, nbGroup_min=100, threshold_potential=0.1, distance_ref=NULL,
         multiV=FALSE, digit.pmin=7, epsilon=10^(-3), epsilon_corrSpat=epsilon*10, iter_max=100,
         trace=TRUE, trace_iter=FALSE, trace_radius=FALSE, export.predicteur=FALSE)
}
\description{
Estimate a finite mixture model with glm-response model, possible concomittant model, two-level spatial regularization and Growing Region step. 
}
\arguments{
  \item{M}{The number of responses variable. \code{integer}.}
  \item{G}{The number of groups to be considered for the mixture. \code{integer}.}
  \item{data}{A dataframe containing the data. \code{data.frame[n,p]}.}
  \item{Id}{The identifier of the images to segment. \code{character vector[n]}. Default is \code{NULL} indicating a unique identifier for all voxels.}
  \item{Wweight_SR}{The local neighborhood matrix. \code{dgCMatrix[n,n]}. Should be normalized by row (i.e. \code{rowSums(Wweight_SR)=1}). Only used if \code{test.GR} or \code{test.ICM} is \code{true}.}
  \item{Wdist_LR}{The regional neighborhood matrix. \code{dgCMatrix[n,n]}. Should contain the distances between the observations (\code{0} indicating infinite distance). Only used if both \code{test.ICM} and \code{test.ICMregional} are \code{true}.}
  
  \item{formula_reg}{A list of formula corresponding to the reponse model of each group for each response variable. \code{list[[M]][[G]]}.}
  \item{offset_reg}{A list of offset corresponding to each reponse model. \code{list[[M]][n,G]}. Default is \code{NULL} indicating no offset.}
  \item{family}{A list of families corresponding to the response model of each group for each response variable. \code{list[[M]][[G]]}.}
  \item{prior_theta}{Initialisation of the means of each group means by random sampling (\code{NULL}), kmeans (\code{"kmeans"}) or user defined values (\code{list[[M]][G]}).}
  \item{prior_sigma}{Initialisation of the standard deviation of each group by dividing the total variance by \code{G} (\code{NULL}), by the standard deviation of the kmeans groups (\code{"kmeans"}) or by user defined values (\code{list[[M]][G]}).}
  \item{posterior_proba}{Initialisation of the posterior membership probabilities to user defined values (\code{matrix[n,g]}). Ignored if \code{NULL}.}
 
  \item{formula_group}{A formula for the concomitant model (\code{formula}). Else (\code{NULL}) indicates no concomittant model. }
  \item{prior_prevalence}{Should a prior based on the prevalence of each group be used ? \code{logical}.}
  \item{prior_proba}{Initialisation of the prior membership probabilities to user defined values (\code{matrix[n,g]}) or to uninformative values (\code{NULL}).}
 
  \item{test.GR}{Should Growing Region algorithm be performed ? \code{logical}.}
  \item{sigma_GR}{Maximum variance of the Growing Region. Can be set to a user defined value (\code{numeric}) or estimated by the variance of group G (\code{"auto"}).}
  \item{proba_GRseed}{Seeds with group probability membership below \code{proba_GRseed} are excluded from the initialisation step of the GR algorithm. \code{numeric}. Ignored if \code{NULL}.}
  \item{proba_GR}{Minimum probability membership to group G required for an observation to be a candidate for the GR. \code{numeric}. Ignored if \code{NULL}.}
  \item{seed}{(\code{T/F}) or (\code{1/0}) indicating the seeds to be used in the GR algorithm. \code{logical vector[n]]}. Can alternatively contains the name of a column in data (\code{character}).}
  \item{n.order}{The penalization for the voxels outside the lesion group defined by the growing region algorithm (bandwith of the gaussian kernel). \code{numeric}. Default is \code{3}.}

  \item{test.ICM}{Should local/regional regularization be performed ? \code{logical}.}
  \item{rho_ICM}{Value of the spatial regularisation parameters : can be fixed a priori (\code{numeric} if test.ICMregional=F and \code{numeric vector[2]} if test.ICMregional=T) or estimated (\code{"init"}).}
  \item{G_ICM}{Potential group merging for the regularization step (\code{interger[G]}) where each element must be inferior or equal to G.}
  \item{prior_prevalenceICM}{Should a prior based on the prevalence of each group be used for estimation of the regularization paramters ? \code{logical}. Default is \code{FALSE}.}
  \item{rho_max}{Maximum possible rho value (\code{numeric}), minimum is 0. Only used if \code{rho_ICM="init"}.}
  \item{update_rho}{Should rho be re-estimated using the posterior probabilites at each step ? \code{logical}.}
 
  \item{test.ICMregional}{Should regional regularization be performed ? \code{logical}. \code{test.ICM} must be also T to be active.}
  \item{coords}{Coordinates (\code{matrix[n,.]}) or the name of columns in data (\code{character vector[.]}) giving the observation coordinates.}
  \item{nbGroup_min}{The minimum group size of the spatial groups required for performing regional regularization. \code{integer}. Default is \code{100}.}
  \item{threshold_potential}{The minimum value of the posterior probability for group G for being considered as lesioned when forming the spatial groups. \code{numeric}.}
  \item{distance_ref}{The distance defining the several neighborhood orders relatively to \code{Wdist_LR}. \code{numeric vector[.]}.}
  \item{multiV}{Should the regional potential range be specific to each spatial group ? \code{logical}. Default is \code{FALSE}.}

  \item{digit.pmin}{Below \code{10^{-digit.pmin}} the posterior probability is set to 0. \code{integer}}
  \item{epsilon}{Convergence occurs when the relative variation of the log-likelihood or of the parameter values between two iterations is below \code{epsilon}. \code{numeric}}
  \item{epsilon_corrSpat}{Spatial regularization begins when the relative variation of the log-likelihood or of the parameter values between two iterations is below \code{epsilon_corrSpat}. \code{numeric}.}
  \item{iter_max}{Maximum number of EM iterations. \code{integer}. Default is \code{100}.}

  \item{trace}{Should initialisation and final EM estimations be displayed ? \code{logical}. Default is \code{TRUE}.}
  \item{trace_iter}{Should estimations be displayed at each EM iteration ? \code{logical}. Default is \code{FALSE}.}
  \item{trace_radius}{Should the radius of the spatial groups be displayed at each EM iteration ? \code{logical}. Default is \code{FALSE}.}
  \item{export.predicteur}{Should the fitted glm models be exported or only their summary ? \code{logical}. Default is \code{FALSE}.}
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
  \item{radius}{if any, a \code{vector} indicating the radius of the spatial groups of group G.}
  \item{group}{a \code{list} indicating for each group the most probable index of observations that belong to the group.}
  \item{family}{a \code{list} containing the family used in the response models.}
  \item{critere_cv}{a \code{list} containing the convergence criteria.}
  \item{cv}{a \code{logical} indicating whether convergence occured.}
  \item{tps}{a \code{list} indicating the time at which the algorithm started and stopped.}  
}
\details{
If a formula \code{formula_group} must not contains any response variable (i.e. must be of the form  \code{~ X1 + X2}).
If \code{formula_group} is \code{NULL}, \code{prior_prevalence} is useless.
If \code{posterior_proba} is specified, other initialisation variables (\code{prior_theta}, \code{prior_sigma}, \code{prior_proba}) are useless.

Growing region algorithm is performed considering only group \code{G}. Pixels that are not included in the result of the growing region algorithm have their probability of belonging to group \code{G} lowered using a gaussian kernel.

\code{G_ICM} indicates the groups used in the regularization step. Default is \code{1:G} which indicates no group merging. For \code{G=4}, \code{G_ICM = c(1,2,2,3)} indicates to merge groups 2 and 3 (i.e. to sum the associated probabilites) in the regularization step.

This function require the following packages : \code{Matrix}, \code{nnet} and \code{Rcpp}, \code{RcppArmadillo} in case of spatial regularization.}
\seealso{
\code{\link{launcher_fMM}} a more user frendly interface to \code{fMMalgo}, \cr
\code{\link{plotCv_fMM}} to plot the convergence criteria, \cr
\code{\link{plotProba_fMM}} to display the prior or the posterior distirbution, \cr
\code{\link{plotY_fMM}} to display the fitted mixture distribution, \cr
\code{\link{calcW_fMM}} to compute the neighborhood matrix, \cr
\code{\link{rhoMF}} to estimate the regularization parameters on external data.
}
\examples{
# see launcher_fMM
}


