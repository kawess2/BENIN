\name{selectFieldDim}
\title{Extract the data dimension}
\alias{selectFieldDim}
\alias{selectFieldDim,Carto3D-method}
\alias{selectFieldDim,MRIaggr-method}
\description{
Extract the dimensions of the spatial field from a \code{\linkS4class{Carto3D}} or from a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{selectFieldDim}{Carto3D}(object)

\S4method{selectFieldDim}{MRIaggr}(object)
}
\arguments{
  \item{object}{an \code{object} of class \code{\linkS4class{Carto3D}} or \code{\linkS4class{MRIaggr}}. REQUIRED.}
}
\value{
  A \emph{data.frame} with one line and three columns named \code{"i","j","k"} indicating the corresponding dimension.
}
\examples{
#### 1- Carto3D method ####
## load NIFTI files and convert them to Carto3D
path.Pat1 <- system.file("nifti", package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file.path(path.Pat1, "TTP_t0"), format = "nifti")
Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0, identifier = "Pat1", param = "TTP_t0")

## selection
selectFieldDim(Carto3D.Pat1_TTP_t0)

#### 2- MRIaggr method ####
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## selection
selectFieldDim(MRIaggr.Pat1_red)
}
\concept{select.}
\keyword{methods}

