\name{initMask}
\title{Initialization of the slice numbers}
\alias{initMask}
\alias{initMask,MRIaggr-method}
\description{
Check and initialize the \code{mask} argument of a \code{\linkS4class{MRIaggr}} method. For internal use.
}
\usage{
\S4method{initMask}{MRIaggr}(object, mask, test = TRUE, init = TRUE,
         format = "data.frame", arg_name = "mask", long_name = "mask", method)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{Carto3D}} or \code{\linkS4class{MRIaggr}}.}
  \item{mask}{the binary contrast parameters to check or initialize. \emph{character vector}.}
  \item{test}{should the mask be checked ? \emph{logical}.  }
  \item{init}{should the mask be initialized ? \emph{logical}. See the details section.}
  \item{format}{the format of the output. Any of \code{"vector"}, \code{"matrix"} and \code{"data.frame"}.}  
  \item{arg_name}{a short name for the error message. \emph{character}.}
  \item{long_name}{the complete name for the error message. \emph{character}.}
  \item{method}{the name of the function that called the initializer. \emph{character}.}
}
\details{
ARGUMENTS : \cr
Setting \code{init} to \code{TRUE} leads to apply \code{numeric2logical} to each mask.
}
\concept{init.}
\keyword{methods,internal}
