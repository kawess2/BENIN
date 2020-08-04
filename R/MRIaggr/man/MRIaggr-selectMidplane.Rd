\name{selectMidplane}
\title{Extract the position of the mid-sagittal plan}
\alias{selectMidplane}
\alias{selectMidplane,MRIaggr-method}
\description{
Extract the position of the mid-sagittal plan from a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{selectMidplane}{MRIaggr}(object)
}
\arguments{
  \item{object}{an \code{object} of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
}
\value{
  A two column \emph{matrix}.
}
\examples{
data("MRIaggr.Pat1_red", package = "MRIaggr")

## selection
selectMidplane(MRIaggr.Pat1_red)
}        
\concept{select.}                                    
\keyword{methods}

