\name{selectCoords}
\title{Extract spatial coordinates}
\alias{selectCoords}
\alias{selectCoords,Carto3D-method}
\alias{selectCoords,MRIaggr-method}
\description{
Extract the coordinates from a \code{\linkS4class{Carto3D}} or from a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{selectCoords}{Carto3D}(object, coords=c("i","j","k"), num = NULL, format = "data.frame")

\S4method{selectCoords}{MRIaggr}(object, coords = c("i","j","k"), spatial_res = c(1,1,1),
         num = NULL, hemisphere = "both", subset = NULL, 
		 slice_var = "k", format = "data.frame")
}
\arguments{
  \item{object}{an \code{object} of class \code{\linkS4class{Carto3D}} or \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{coords}{the coordinates that sould be extracted. Any of \code{"i"} \code{"j"} \code{"k"} or \code{"index"}.}
  \item{spatial_res}{a dilatation factor for the coordinates. \emph{positive numeric vector of size 3}.}
  \item{num}{the slices to extract. \emph{numeric vector} or \code{NULL}.}
  \item{hemisphere}{the hemisphere to extract. \emph{character}.}
  \item{subset}{the subset of observations to extract. \emph{positive integer vector} or \code{NULL} leading to use all observations}
  \item{slice_var}{the type of slice to extract. \code{"i"} for sagittal, \code{"j"} for coronal and \code{"k"} for transverse. \emph{character}.}  
  \item{format}{the format of the output. Can be \code{"vector"}, \code{"matrix"} or \code{"data.frame"}.}
}
\details{
ARGUMENTS: \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

Information about the \code{hemisphere} argument can be found in the details section of \code{\link{selectContrast}}.

FUNCTION: \cr
Each of the \code{num}, \code{hemisphere} and \code{subset} argument define a subset of the total set of observations. 
It is the intersection of all these three subsets that is extracted.
}
\value{
  A \code{"vector"}, \emph{matrix} or a \emph{data.frame}. In the latter two cases, each row corresponds to a voxel and each column to a coordinate.
}
\seealso{
 \code{\link{calcHemisphere}} to identify the hemispheres. \cr
 \code{\link{allocHemisphere<-}} and \code{\link{allocContrast<-}} to allocate hemispheres obtained from an external source.
}
\examples{
#### 1- Carto3D method ####
## load nifti files and convert them to Carto3D
path.Pat1 <- system.file("nifti", package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file.path(path.Pat1, "TTP_t0"), format = "nifti")
Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0, identifier = "Pat1", param = "TTP_t0")

## selection all coordinates
coords1 <- selectCoords(Carto3D.Pat1_TTP_t0)
dim(coords1)

## selection coordinates i and j from slices 1 to 3
coords2 <- selectCoords(Carto3D.Pat1_TTP_t0, num = 1:3, coords = c("i","j"))
dim(coords2)

#### 2- MRIaggr method ####
## load a MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## select all coordinates for all observations
coords <- selectCoords(MRIaggr.Pat1_red)
dim(coords)
head(coords)

## select coordinate i for slices 1 and 3
coords <- selectCoords(MRIaggr.Pat1_red, coords = "i", num = c(1,3))
dim(coords)
head(coords)

## select coordinate i for observations in the hemishere containing the lesion
coords <- selectCoords(MRIaggr.Pat1_red, hemisphere = "lesion", num = c(1,3))
dim(coords)
head(coords)

## select coordinate i for observations in the right hemisphere
coords <- selectCoords(MRIaggr.Pat1_red, hemisphere = "right", num = c(1,3))
dim(coords)
head(coords)

## select all coordinates and rescale them
coords <- selectCoords(MRIaggr.Pat1_red, spatial_res = c(1.875,1.875,6))
dim(coords)
head(coords)

## select coordinate i and j and return a matrix
coords <- selectCoords(MRIaggr.Pat1_red, format = "matrix")
is(coords)
head(coords)

}
\concept{select.}
\keyword{methods}

