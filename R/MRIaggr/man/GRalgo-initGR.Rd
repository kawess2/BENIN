\name{initGR}
\title{Growing Region initialization}
\alias{initGR}
\usage{
initGR(contrast, W, seed, sigma_max, range, range.seed, breaks, iter_max, sd.robust,
         keep.lower, keep.upper, history.sigma, history.step, history.front,
         rescale, verbose, method = "initGR")
}
\description{
Initialize arguments of the Growing Region algorithm. For internal use.
}
\arguments{
  \item{contrast}{the contrast value of each observation. \emph{numeric vector}.}
  \item{W}{the neighbourhood matrix. \emph{dgCMatrix}.}
  \item{seed}{the index of the initial seeds or a binary indicator of the initial seeds. \emph{positive integer vector} or \emph{logical vector}.}
  \item{sigma_max}{the maximum admissible value for the variability of the group contrast. \emph{positive numeric}.}
  \item{range}{the range of acceptable contrast values for the growing region group. \emph{numeric vector of size 2}.}
  \item{range.seed}{the range of acceptable contrast values for the seeds. \emph{numeric vector of size 2}.}
 
  \item{breaks}{the break points or the number of break points to use to categorize the contrast distribution. \emph{numeric vector} or \code{postive integer}.}
  \item{iter_max}{the maximum number of iterations for the expansion of the growing region. \code{postive integer}.}  
  \item{sd.robust}{should the median absolute deviation be used to estimte the variability of the group contrast, or the standard deviation ? \emph{logical}.}
  \item{keep.lower}{should removing observations with high intensity of the region be forbidden ? \emph{logical}.}  
  \item{keep.upper}{should removing observations with low intensity of the region be forbidden ? \emph{logical}.}  
 
  \item{history.sigma}{should the values of sigma be recorded ? \emph{logical}.}
  \item{history.step}{should the number of observations included in the growing region set be recorded ? \emph{logical}.}
  \item{history.front}{should the propagation front of the GR set be recorded ? \emph{logical}.}

  \item{rescale}{should the contrast be scaled ? \emph{logical}.}  
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}  
  \item{method}{the name of the function that called the initializer. \emph{character}.}
}
\concept{init.}
\keyword{functions,internal}
