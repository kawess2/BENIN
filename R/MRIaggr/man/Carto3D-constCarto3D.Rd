\name{constCarto3D}
\title{Array constructor for Carto3D objects}
\alias{constCarto3D}
\usage{
constCarto3D(array, identifier, param, default_value = NULL,            
         pos_default_value = c(1,1,1), voxelDim = NULL, rm.array = FALSE)
}
\description{
Creates a \code{\link{Carto3D}} object from an array.
}
\arguments{
  \item{array}{the values of the contrast parameter. \code{array}. REQUIRED.}
  \item{identifier}{the identifier of the patient from which the data originated. \emph{character}. REQUIRED.}
  \item{param}{the contrast parameter. \emph{character}. REQUIRED.}
  \item{default_value}{the reference values of the contrast parameters (e.g. the background value). \emph{character} or \code{NULL} leading to search the reference value in \code{array[pos_default_value]}.}
  \item{pos_default_value}{the coordinates of the observation that contains the reference value. \emph{numeric vector} with length the number of dimension of array.}
  \item{voxelDim}{the voxel size with its unit of measure. A four columns \emph{data.frame} with names \code{"i"}, \code{"j"}, \code{"k"} and \code{"unit"}.}
  \item{rm.array}{should the object on which \emph{array} argument points be removed form the global environment ? \emph{logical}.}
}
\details{
ARGUMENTS: \cr
\code{pos_default_value} is active only if \code{default_value} is set to \code{NULL}.
}
\value{
  An object of class \code{\linkS4class{Carto3D}}.
}
\examples{
## load NIFTI files
path.Pat1 <- system.file(file.path("nifti"), package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file.path(path.Pat1, "TTP_t0"), format = "nifti")

## convert them to Carto3D
Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0, identifier = "Pat1", param = "TTP_t0")
class(Carto3D.Pat1_TTP_t0)
}
\concept{const.}
\keyword{functions}

