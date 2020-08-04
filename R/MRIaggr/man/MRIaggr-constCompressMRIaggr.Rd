\name{constCompressMRIaggr}
\title{Compress a MRIaggr object}
\alias{constCompressMRIaggr}
\alias{constCompressMRIaggr,MRIaggr-method}
\description{
Generate from an existing \code{\linkS4class{MRIaggr}} object a new one with lower spatial resolution.
}
\usage{
\S4method{constCompressMRIaggr}{MRIaggr}(object, compression.factor, param = NULL,
         mask = NULL, threshold = 0.49, verbose = optionsMRIaggr("verbose"))
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{compression.factor}{the compression factor. \emph{postive integer}. REQUIRED.}
  \item{param}{the contrast parameters to load in the new MRIaggr object. \emph{character vector} or \code{NULL}.}
  \item{mask}{the binary contrast parameter(s). \emph{character vector}.}
  \item{threshold}{the value above which the local mean of the binary parameters is assigned to 1 (and otherwise to 0). \emph{numeric between 0 and 1}.}
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
}
\value{
  a \code{\linkS4class{MRIaggr}} object.
}
\details{
ARGUMENTS: \cr
Information about the \code{param} argument can be found in the details section of \code{\link{initParameter}}.

FUNCTION: \cr
The fonction uses a local mean to compress the initial parameters maps to a lower resolution. 
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")
 
## compress the MRIaggr object
MRIaggr.compressed <- constCompressMRIaggr(MRIaggr.Pat1_red, compression.factor = 2,
                                           param = c("DWI_t0","T2_FLAIR_t2","MASK_T2_FLAIR_t2"),
                                           mask = "MASK_T2_FLAIR_t2") 

## display
par(mfrow = c(2,4), mar = rep(1.75,4), mgp = c(2,0.75,0))
multiplot(MRIaggr.Pat1_red, param = "DWI_t0", window = NULL, breaks = seq(0,350,1),
             midplane = TRUE, main = "before - slice ")
multiplot(MRIaggr.compressed, param = "DWI_t0", window = NULL, breaks = seq(0,350,1),
             midplane = TRUE, main = "after - slice ")

multiplot(MRIaggr.Pat1_red, param = "MASK_T2_FLAIR_t2", main = "before - slice ")
multiplot(MRIaggr.compressed, param = "MASK_T2_FLAIR_t2", main = "after - slice ")
}
\concept{const.}
\keyword{methods}

