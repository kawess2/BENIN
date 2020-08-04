\name{initDisplayWindow}
\title{Device management}
\alias{initDisplayWindow}
\description{
Display the device on which the plot will be displayed. For internal use.
}
\usage{
initDisplayWindow(window, filename, path, width, height, scale, res,
         mfrow, bg, pty, mar, mgp, n.contrast = 1)
}
\arguments{
  \item{window}{the type of device on which the plot will be displayed. \emph{logical}, \code{NULL} or \emph{character}.}
  \item{filename}{the name of the file used to export the plot. \emph{character}.}
  \item{path}{the directory where the plot file will be created. \emph{character}.}
  \item{width}{the width of the device used to export the plot in inches. \emph{postive numeric}.}
  \item{height}{the height of the device used to export the plot. \emph{postive numeric}.}
  \item{scale}{the scaling factor to convert \code{height} and \code{height} to standard unit. \emph{numeric}.}
  \item{res}{the nominal resolution in ppi which will be recorded in the bitmap file. \emph{positive integer}.}

  \item{mfrow}{the division of the device in plot region. \emph{numeric vector of size 2}.}
  \item{bg}{the color used for the background. \emph{character}.}
  \item{pty}{the type of plot region to be used. Can be \code{"s"} or \code{"m"}.}
  \item{mar}{the number of margin lines to be specified on the four sides of the plot. \emph{positive numeric vector of size 4}.}
  \item{mgp}{the margin line for the axis title, axis labels and axis line. \emph{positive numeric vector of size 3}.}
  \item{n.contrast}{the number of contrast parameters. \emph{postive integer}.}
}
\details{
ARGUMENTS: \cr
Information about the \code{window}, \code{filename}, \code{width}, \code{height}, \code{path}, \code{unit} and \code{res} arguments can be found in the details section of \code{\link{initWindow}}.

Information about the \code{bg}, \code{pty}, \code{mar} and \code{mgp} arguments can be found in \code{\link{par}}.
}
\concept{init.}
\keyword{function,internal}
