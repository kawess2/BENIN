\name{selectIdentifier}
\title{Extract the identifier}
\alias{selectIdentifier}
\alias{selectIdentifier,Carto3D-method}
\alias{selectIdentifier,MRIaggr-method}
\description{
Extract the patient identifier from a \code{\linkS4class{Carto3D}} or from a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{selectIdentifier}{Carto3D}(object)

\S4method{selectIdentifier}{MRIaggr}(object)
}
\arguments{
  \item{object}{an \code{object} of class \code{\linkS4class{Carto3D}} or \code{\linkS4class{MRIaggr}}. REQUIRED.}
}
\value{
  A \emph{character}.
}
\examples{
#### 1- Carto3D method ####
## load nifti files and convert them to Carto3D
path.Pat1 <- system.file("nifti", package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file.path(path.Pat1, "TTP_t0"), format = "nifti")
Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0, identifier = "Pat1", param = "TTP_t0")

## selection
selectIdentifier(Carto3D.Pat1_TTP_t0)

#### 2- MRIaggr method ####
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## selection
selectIdentifier(MRIaggr.Pat1_red)
}        
\concept{select.}                                    
\keyword{methods}

