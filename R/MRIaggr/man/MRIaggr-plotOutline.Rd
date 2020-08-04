\name{pointsOutline}
\title{Compute the outline of a spatial group}
\alias{pointsOutline}
\usage{
pointsOutline(coords, array = NULL, filter = "2D_N4")
}
\description{
Return a the outline of a spatial group for a given neighbourhood. For internal use.
}
\arguments{
  \item{coords}{the spatial coordinates of the observations. \emph{data.frame}.}
  \item{array}{alternative specification of the spatial coordinates using an array where the non-NA values indicates the points of interest. \emph{array} or \code{NULL} leading to use the \code{coords} argument.}
  \item{filter}{the type of filter, see \code{\link{calcFilter}} for more details. \emph{character}.}
}

\details{
FUNCTION: \cr
 This function apply a filter to the spatial group and consider as the outline all the observations that do not reach the maximum value (i.e. do not have a complete neighbourhood). 
}
\value{
  A \emph{data.frame} with the spatial coordinates of the outline.
}
\concept{calc.}
\keyword{functions,internal}