\name{selectVoxelDim}
\title{Extract the dimensions of a voxel}
\alias{selectVoxelDim}
\alias{selectVoxelDim,MRIaggr-method}
\description{
Extract the voxel size stored in the \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{selectVoxelDim}{MRIaggr}(object, unit = TRUE)
}
\arguments{
  \item{object}{an \code{object} of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{unit}{should the unit be returned ? \emph{logical}.}
}
\value{
  A \emph{data.frame} with one line and three or four columns containing the voxel size for each dimension and the size unit (if requested).
}
\examples{
data("MRIaggr.Pat1_red", package = "MRIaggr")

## selection
selectVoxelDim(MRIaggr.Pat1_red)
}
\concept{select.}
\keyword{methods}

