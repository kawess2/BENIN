\name{plotDistClass}
\title{Plot the distribution of the contrast parameter}
\alias{plotDistClass}
\alias{plotDistClass,MRIaggr-method}
\description{
Distribution of a contrast parameter according to the tissue class.
}
\usage{
\S4method{plotDistClass}{MRIaggr}(object, param, param.membership, num = NULL,
         bw.adjust = 1, kernel = "gaussian", from = NULL, to = NULL, ylim = NULL,
		 col = 1:6, main = NULL, mgp = c(2,0.5,0), type = "l",
         pch = 20, lwd = 1, x.legend = "topright", y.legend = NULL, cex.legend = 0.8,
         filename = paste(object@identifier, "plotDistClass", sep = "_"), ...)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{param}{the contrast parameter to display. \emph{character}. REQUIRED.}
  \item{param.membership}{the parameters indicating the probabilistic membership to the tissue classes. \emph{character vector}.}
  \item{num}{the slices to use. \emph{numeric vector} or \code{NULL}.}

  \item{bw.adjust}{the smoothing bandwidth to use. \emph{numeric}. See \code{\link{density}} for more details.}
  \item{kernel}{the smoothing kernel to use. \emph{character}. See \code{\link{density}} for more details.}
  \item{from,to}{the left and right-most points of the grid at which the density is to be estimated. \emph{numeric} or \code{NULL} leading to automatic adjustment. See \code{\link{density}} for more details.}
  \item{ylim}{the y limits of the plot. \emph{numeric vector of size 2} or \code{NULL} leading to automatic setting of the y limits.}

   \item{col}{the colors with which the distributions will be displayed. \emph{character vector}.}
  \item{main}{an overall title for the plot. \emph{character}.}
  \item{mgp}{the margin line for the axis title, axis labels and axis line. \emph{positive numeric vector of size 3}.}
  \item{type}{the type of plot to display. \emph{character}. See \code{\link{plot.default}} for more details.}
  \item{pch}{the symbol with which the distribution will be displayed. \emph{positive integer}.}
  \item{lwd}{the line width. \emph{postive numeric}.}
  \item{x.legend}{the x coordinates of the legend. \emph{numeric} or \emph{character}.}
  \item{y.legend}{the y coordinates of the legend. \emph{numeric} or \emph{character}.}
  \item{cex.legend}{the expansion factor of the legend. \emph{positive numeric}.}		

  \item{filename}{the name of the file used to export the plot. \emph{character}.}
  \item{...}{additional arguments to be passed to \code{\link{optionsMRIaggr}} for specifying the graphical parameters.}
}
\details{
ARGUMENTS: \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

Information about the \code{lwd} argument can be found in \code{\link{par}}.

Information about the \code{x.legend}, \code{y.legend}, \code{cex.legend} arguments can be found in \code{\link{legend}} (\code{cex.legend} is the \code{cex} argument of \code{legend}).

Argument(s) \code{...} must correspond to some of the following arguments : \code{height}, \code{hemisphere}, \code{norm_mu}, \code{norm_sigma}, \code{path}, \code{res}, \code{unit} \code{width}, \code{window}.

FUNCTION: \cr
This method relies on the \code{\link{density}} function.
}
\value{
  None.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## display
plotDistClass(MRIaggr.Pat1_red, param = "DWI_t0",
              param.membership = c("MASK_T2_FLAIR_t2"))

# specify the smoothing bandwidth
plotDistClass(MRIaggr.Pat1_red, param = "DWI_t0", bw.adjust = 1,
              param.membership="MASK_T2_FLAIR_t2")

# specify the scale
plotDistClass(MRIaggr.Pat1_red, param = "DWI_t0", bw.adjust = 1,
              from = 200, to = 300, param.membership = "MASK_T2_FLAIR_t2")

# use several classes
\dontrun{	   
plotDistClass(MRIaggr.Pat1_red, param = "TTP_t0", bw.adjust = 2,
              param.membership = c("CSF","WM","GM","MASK_T2_FLAIR_t2"))

plotDistClass(MRIaggr.Pat1_red, param = "DWI_t0", bw.adjust = 2,
              param.membership = c("CSF","WM","GM","MASK_T2_FLAIR_t2"))
}

}
\concept{plot.}
\keyword{methods}

