\name{selectN}
\title{Extract the number of observations}
\alias{selectN}
\alias{selectN,Carto3D-method}
\alias{selectN,MRIaggr-method}
\description{
Extract the number of observations contained in a \code{\linkS4class{Carto3D}} or in a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{selectN}{Carto3D}(object, num = NULL)
		 
\S4method{selectN}{MRIaggr}(object, num = NULL, hemisphere = "both", subset = NULL)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{num}{the slices to consider. \emph{numeric vector} or \code{NULL}.}
  \item{hemisphere}{the hemisphere to consider. \emph{character}.}
  \item{subset}{the subset of observations to consider. \emph{positive integer vector} or \code{NULL} leading to consider all observations.}
}
\details{
ARGUMENTS: \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

Information about the \code{hemisphere} argument can be found in the details section of \code{\link{selectContrast}}.

FUNCTION: \cr
Each of the \code{num}, \code{hemisphere} and \code{subset} argument define a subset of the total set of observations. 
It is the length of the intersection of all these three subsets that is measured.
}
\value{
  An \emph{integer}.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## total number of observations
res <- selectN(MRIaggr.Pat1_red)

## number of observations for the hemisphere that contains the lesion
res <- selectN(MRIaggr.Pat1_red, hemisphere = "lesion")

## number of observations in the first 1000 observations t
## that are in the hemisphere containing the lesion
res <- selectN(MRIaggr.Pat1_red, subset = 1:1000, hemisphere = "lesion")
}
\concept{select.}
\keyword{methods}

