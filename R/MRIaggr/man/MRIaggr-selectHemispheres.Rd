\name{selectHemispheres}
\title{Extract the position of the lesion in each hemisphere}
\alias{selectHemispheres}
\alias{selectHemispheres,MRIaggr-method}
\description{
Extract the position of the lesion in each hemisphere from a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{selectHemispheres}{MRIaggr}(object, hemisphere = "both")
}
\arguments{
  \item{object}{an \code{object} of class \code{\linkS4class{Carto3D}} or \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{hemisphere}{the hemisphere of interest. \emph{character}. See the details section.}
}
\value{
  A \emph{character vector} or a \emph{data.frame}.
}
\details{
Possible values for the \code{hemisphere} argument are:
\itemize{
\item \code{"both"} : indicates the presence of the lesion in both hemispheres.
\item \code{"left"} : indicates the presence of the lesion in the left hemisphere.
\item \code{"right"} : indicates the presence of the lesion in the right hemisphere.
\item \code{"lesion"} : indicates which hemispheres are lesioned. 
\item \code{"contralateral"} : indicates which hemispheres are intact.
}
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## selection
selectHemispheres(MRIaggr.Pat1_red)
selectHemispheres(MRIaggr.Pat1_red, hemisphere = "right")
selectHemispheres(MRIaggr.Pat1_red, hemisphere = "left")
selectHemispheres(MRIaggr.Pat1_red, hemisphere = "lesion")
selectHemispheres(MRIaggr.Pat1_red, hemisphere = "contralateral")
}        
\concept{select.}                                    
\keyword{methods}

