\name{calcTissueType}
\title{Probabilistic tissue type segmentation}
\alias{calcTissueType}
\alias{calcTissueType,MRIaggr-method}
\description{
Perform a probabilistic segmentation of the voxel in Cerebro Spinal Fluid, White Matter and Grey Matter classes.
}
\usage{
\S4method{calcTissueType}{MRIaggr}(object, param, niter = 100, nnei = 6,
         beta = if(sub == TRUE){0.3}else{0.7}, sub = TRUE, digit = 0, verbose = TRUE,
         name_newparam = c("CSF", "GM", "WM"), update.object = FALSE, overwrite = FALSE)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{param}{the contrast parameter that should be used to distinguish the WM, the GM and the CSF. \emph{character}. REQUIRED.}
  \item{niter}{the number of iterations used by \code{mritc.bayes}. \emph{positive integer}.}
  \item{nnei}{the number of neighbors. \emph{positive integer}.}
  \item{beta}{the parameter 'inverse temperature' of the Potts model. \emph{numeric}.}
  \item{sub}{if \code{TRUE}, use the higher resolution model; otherwise, use the whole voxel method. \emph{logical}.}
  \item{digit}{the number of decimal places to use for the initialization. \emph{positive integer}.}

  \item{verbose}{indicate the level of output as the algorithm runs. \emph{logical}.}
  \item{name_newparam}{the name of the new paramaters containing the probabilistic segmentation. \emph{character vector of size 3}.}  
  \item{update.object}{should the resulting tissue types be stored in \code{object} ? \emph{logical}.}
  \item{overwrite}{if tissue types are already stored in \code{object@data}, can they be overwritten ? \emph{logical}.}  
}
\details{

This function requires to have installed the \emph{mritc} package to work.

ARGUMENTS: \cr
Information about the \code{nnei} and \code{sub} arguments can be found in \code{makeMRIspatial}. \cr
Information about the \code{niter} and \code{beta} arguments can be found in \code{mritc}.
  
FUNCTION: \cr
This function uses the \code{mritc.bayes} function of the \emph{mritc} package to compute the probabilistic segmentation.
T1 sequence is the recommanded sequence to identify the various tissue types but T2 gradient echo may also be used.

The initialization function \code{initOtsu} was found much more slower when the contrast parameter values have a large number of digit. The \code{digit} argument enable to round the contrast values ONLY for the computation of the initialization values.
}
\value{
  An \emph{list} containing :
  \itemize{
  \item\code{[[prob]]} : the group membership of each voxel for each class. \emph{matrix}.
  \item\code{[[mu]]} : the mean value of each class. \emph{numeric vector}.
  \item\code{[[sigma]]} : the standard deviation of each class. \emph{numeric vector}.
  }
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

\dontrun{
## perform segmentation  (call mritc)
calcTissueType(MRIaggr.Pat1_red, param = "T1_t0", update.object = TRUE, overwrite = TRUE)

## display
multiplot(MRIaggr.Pat1_red, num = 1,
             param = c("CSF","WM","GM"), legend = FALSE,
             palette = "rgb")
}
}
\concept{calc.}
\keyword{methods}

