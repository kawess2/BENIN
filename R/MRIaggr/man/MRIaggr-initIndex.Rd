\name{initIndex}
\title{Index initialization}
\alias{initIndex}
\description{
Check and initialize index. For internal use.
}
\usage{
initIndex(object, index, num, hemisphere = "both", numeric2logical, indexNum = NULL,
         outline.default, cex.default, pch.default, col.default, filter.default, method)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}} or \code{NULL} .}
  \item{index}{the coordinates of additionnal points to display. \emph{data.frame} or \emph{list}.}
  \item{num}{the slices to display. \emph{numeric vector} or \code{NULL}.}
  \item{hemisphere}{the hemisphere to display. \emph{character}.}
  \item{numeric2logical}{if a parameter is specified for the index arguments, should it be converted to logical ? \emph{logical}.}
  \item{indexNum}{the number associated to the index (for display). \emph{numeric}.}
  \item{outline.default}{the default display for the spatial group (edge or entire group). \emph{logical}.}
  \item{cex.default}{the default expansion factor used to plot the observations. \emph{numeric}.}
  \item{pch.default}{the default label used to plot the observations. \emph{numeric}.}
  \item{col.default}{the default color used to plot the observations. \emph{character vector}.}
  \item{filter.default}{the default filter used to define the neighbourhood. \emph{character}.} 
  \item{method}{the name of the function that called the initializer. \emph{character}.}
}
\details{
ARGUMENTS: \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

Information about the \code{hemisphere} argument can be found in the details section of \code{\link{selectContrast}}.

Possible values for \code{index} are:
\itemize{
\item \code{NULL} : the argument is skiped, no additional points are displayed.
\item \emph{character} : the name of a logical parameter contained in the object (if called by an \code{MRIaggr} method)
\item \emph{data.frame} with 3 columns named \code{"i"},\code{"j"},\code{"k"} containing the coordinates of the points to display.
\item \emph{list} with an element named \code{"data"} containing the coordinates of the points to display (see previous point) or the name of a logical parameter (if called by an \code{MRIaggr} method). \cr
It can also contain  \code{"cex"}, \code{"pch"} and \code{"col"} component to specify the how the points should be displayed. \cr
To only display the outline of the spatial group instead of the spatial group itself, the list must contains an element named \code{"outline"} that have for value \code{TRUE}.  \cr
The type of neighbourhood used to determine the outline can be specified with a \code{"filter"} element in the list. The default neighbourhood is \code{"2D_N4"}.
}
}
\concept{init.}
\keyword{function,internal}
