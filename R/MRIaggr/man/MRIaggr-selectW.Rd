\name{selectW}
\title{Extract a neighbourhood matrix}
\alias{selectW}
\alias{selectW,MRIaggr-method}
\description{
Extract elements in the \code{W} slot of a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{selectW}{MRIaggr}(object, type = "Wmatrix", subset_W = NULL, 
         hemisphere = "both", num = NULL, slice_var = "k", upper = NULL)
}
\arguments{
  \item{object}{an \code{object} of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{type}{the type of the element to select. Can be \code{"Wmatrix"}, \code{"Wblocks"} or \code{"upper"}.}
  \item{subset_W}{the subset of observations to select. \emph{integer vector} or \code{NULL} leading to select all the observations.}
  \item{hemisphere}{the hemisphere to use. \emph{character}.}
  \item{num}{the slices to use. \emph{numeric vector} or \code{NULL}.}
  \item{slice_var}{the type of slice to extract. \code{"i"} for sagittal, \code{"j"} for coronal and \code{"k"} for transverse. \emph{character}.}  
  \item{upper}{the format of the matrix to extract. Can be the full matrix (\code{upper=NULL}) or upper- or lower-triangular matrix (\code{upper=TRUE} or \code{upper=FALSE}).}  
}
\details{
This function requires to have installed the \emph{Matrix} and the \emph{spam} package to work when argument \code{name} is set to \code{"W_euclidean"}.

ARGUMENTS: \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

Information about the \code{hemisphere} argument can be found in the details section of \code{\link{selectContrast}}.

FUNCTION: \cr
If \code{subset_W} is not \code{NULL} then arguments \code{hemisphere} and \code{num} are ignored. \cr
Each of the \code{num}, \code{hemisphere} and \code{subset} argument define a subset of the total set of observations. 
It is the intersection of all these three subsets is extracted.
}
\seealso{ 
 \code{\link{calcW}} to compute the neighboring matrix. \cr
 \code{\link{calcBlockW}} to decompose the neighboring matrix in independant blocks. \cr
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## compute and allocate a neighbourhood matrix
calcW(object=MRIaggr.Pat1_red, range = sqrt(2), update.object = TRUE)

## select the neighbourhood matrix
selectW(MRIaggr.Pat1_red)[1:10,1:10]
selectW(MRIaggr.Pat1_red, upper = TRUE)[1:10,1:10]
# selectW(MRIaggr.Pat1_red, upper = FALSE)[1:10,1:10] # will not work
}
\concept{select.}
\keyword{methods}

