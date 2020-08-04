\name{simulPotts_cpp}
\title{Potts model simulation}
\alias{simulPotts_cpp}
\description{
C++ function that simulates an isotropic Potts model by Gibbs sampling. Experimental.
}
\usage{
simulPotts_cpp(W_SR, W_LR, sample, coords, site_order, rho, 
         distance_ref, iter_max, cv_criterion, regional, verbose)
}
\arguments{
  \item{W_SR}{the local neighbourhood matrix. \emph{dgCMatrix}. Should be normalized by row (i.e. \code{rowSums(W_SR)=1}).}
  \item{W_LR}{the regional neighbourhood matrix. \emph{dgCMatrix}. Should contain the categroy of distance.}
  \item{sample}{the initial observations (by rows) for the various groups (by columns). \emph{matrix}.}  
  \item{coords}{the voxels coordinates. \emph{matrix}.}
  \item{site_order}{a specific order to go all over the sites. An \emph{integer vector} or a negative number to indicate no specific order.}  
  \item{rho}{a \emph{numeric} corresponding to the Potts model parameter (i.e. intensity of the neighbourhood correlations).}
  \item{distance_ref}{the distance defining the several neighbourhood orders relatively to \code{Wdist_LR}. \emph{numeric vector}.}
  \item{iter_max}{the maximum number of iteration of the Gibbs sampler. \emph{integer}.}
  \item{cv_criterion}{if the maximum absolute difference of the probabilistic membership between two iterations is below \code{critere} then the Gibbs sampler is stopped. \emph{numeric}.}
  \item{regional}{should the regional potential be computed ?. \emph{logical}.}
  \item{verbose}{should the convergence criterion be printed at each step ? \emph{logical}.}
}
\concept{cpp.}
\keyword{functions,internal}

                        
                        
                      