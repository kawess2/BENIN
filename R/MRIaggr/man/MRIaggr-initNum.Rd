\name{initNum}
\title{Initialization of the slice numbers}
\alias{initNum}
\alias{initNum,Carto3D-method}
\alias{initNum,MRIaggr-method}
\description{
Check and initialize the \code{num} argument of a \code{\linkS4class{Carto3D}} or a \code{\linkS4class{MRIaggr}} method. For internal use.
}
\usage{
\S4method{initNum}{Carto3D}(object, num, test = TRUE, init = TRUE, method)

\S4method{initNum}{MRIaggr}(object, num, test = TRUE, init = TRUE, slice_var = "k", method)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{Carto3D}} or \code{\linkS4class{MRIaggr}}.}
  \item{num}{the slice numbers to check or initialize. \emph{numeric vector} or \code{NULL}. See the details section.}
  \item{test}{should the slice numbers be checked ? \emph{logical}.  }
  \item{init}{should the slice numbers be initialized if \code{num} equals \code{NULL} ? \emph{logical}.  }
  \item{slice_var}{the type of slice to extract. \code{"i"} for sagittal, \code{"j"} for coronal and \code{"k"} for transverse. \emph{character}.}  
  \item{method}{the name of the function that called the initializer. \emph{character}.}
}
\details{
ARGUMENTS : \cr
Setting \code{num} to \code{NULL} leads to load all available slices.
}
\concept{init.}
\keyword{methods,internal}
