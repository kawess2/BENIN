\name{initCol}
\title{Color initialization}
\alias{initCol}
\description{
Check and initialize display arguments. For internal use.
}
\usage{
initCol(contrast, coords, param = NULL, pch, col, palette,
        breaks, legend, type.breaks, method)
}
\arguments{
  \item{contrast}{the contrast value of each observation. \emph{matrix}.}
  \item{coords}{the spatial coordinates of the observations. \emph{data.frame}.}
  \item{param}{the contrast parameter to display. \emph{character}.}
  \item{pch}{the symbol with which the observations will be displayed. \emph{positive integer}.}
  \item{col}{the color with which the observations will be displayed. \emph{character vector}.}
  \item{palette}{the colors or the palette to use when associating colors to intensities. \emph{character vector} or \emph{character}.}
  \item{breaks}{the break points or the number of breakpoints to use to generate the color intervals. \emph{numeric vector} or \emph{postive integer}.}
  \item{legend}{how should the legend be displayed ? \emph{logical}, \code{NULL} or \code{"only"}.}
  \item{type.breaks}{should the break points be equally space according the range of data values (\code{"range"}), centered (\code{"range_center"}) or correspond to the quantile values (\code{"quantile"}) ?}
  \item{method}{the name of the function that called the initializer. \emph{character}.}
}
\details{
FUNCTION: \cr
\code{breaks} and \code{palette} are active only if \code{col} is \code{NULL}. 
}
\concept{init.}
\keyword{function,internal}
