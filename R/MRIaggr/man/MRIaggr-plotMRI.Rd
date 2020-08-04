\name{plotMRI}
\title{Diplay a contrast parameter by coordinates}
\alias{plotMRI}
\description{
Display a single slice of a contrast parameter. For internal use.
}
\usage{
plotMRI(contrast, coords, breaks, palette, col, asp,
         xlim, ylim, pch, cex, axes, col.NA, pch.NA, xlab, ylab, main, cex.main)
}
\arguments{
  \item{contrast}{the intensities to display. \emph{numeric vector}.}
  \item{coords}{the spatial coordinates of the observations. \emph{data.frame}.}
  \item{breaks}{the break points to use to generate the color intervals. \emph{numeric vector}.}
  \item{palette}{the colors or the palette to use when associating colors to observation intensities. \emph{character vector}.}
  \item{col}{the colors with which the observations will be displayed. \emph{character vector}.}
  \item{asp}{the aspect ratio y/x. \emph{numeric}.}

  \item{xlim}{the x limits of the plot. \emph{numeric vector of size 2}.}
  \item{ylim}{the y limits of the plot. \emph{numeric vector of size 2}.}
  \item{pch}{the symbol with which the observations will be displayed. \emph{positive integer}.}
  \item{cex}{the expansion factor for the observation labels. \emph{positive numeric}.}
  \item{axes}{should the axes be plotted ? \emph{logical}.}
  \item{col.NA}{the color to use to plot the NAs. \emph{character}.}
  \item{pch.NA}{the label to use to plot the NAs. \emph{postive integer}.}
  \item{xlab}{a title for the x axis. \emph{character}.}
  \item{ylab}{a title for the y axis. \emph{character}.}
  \item{main}{an overall title for the plot. \emph{character}.}
  \item{cex.main}{the expansion factor for the main title. \emph{numeric}.}
}
\details{
ARGUMENTS: \cr
Information about the \code{palette} argument can be found in the details section of \code{\link{initCol}}.
}
\concept{init.}
\keyword{function,internal}

