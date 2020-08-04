\name{plotTableLesion}
\title{Lesion volume displayed by slices}
\alias{plotTableLesion}
\alias{plotTableLesion,MRIaggr-method}
\description{
Display the lesion volume by slices using the \code{table_lesion} slot.
}
\usage{
\S4method{plotTableLesion}{MRIaggr}(object, mask, num = NULL, type = "matplot",
         col = 1:5, lty = 1:5, lwd = 1, mgp = c(2,0.5,0), mar = rep(3,4),
         main = paste("lesion - ", object@identifier, sep = ""),
         cex.legend = 1, cex.main = 1, cex.axis = 1, cex.lab = 1,
         filename = paste(object@identifier, "plotTableLesion", sep = "_"), ...)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{mask}{the binary contrast parameter indicating the lesion. \emph{character}. REQUIRED.}
  \item{num}{the slices to display. \emph{numeric vector} or \code{NULL}.}
  \item{type}{the type of plot to display. Can be \code{"matplot"} or \code{"evolution"}.}

  \item{col}{the colors with which the volumes will be displayed. \emph{character vector} or \emph{numeric vector}.}
  \item{lty}{the line type used to represent the volume. \emph{numeric vector}.}
  \item{lwd}{the line width. \emph{postive numeric}.}
  \item{main}{an overall title for the plot. \emph{character}.}
  \item{mgp}{the margin line for the axis title, axis labels and axis line. \emph{positive numeric vector of size 3}.}
  \item{mar}{the number of margin lines to be specified on the four sides of the plot. \emph{positive numeric vector of size 4}.}
  \item{cex.main}{the expansion factor for the main title. \emph{numeric}.}
  \item{cex.legend}{the expansion factor of the legend. \emph{positive numeric}.}		
  \item{cex.axis}{the magnification to be used for axis annotation relative to the current setting of cex. \emph{positive numeric}.}		
  \item{cex.lab}{the magnification to be used for x and y labels relative to the current setting of cex. \emph{positive numeric}.}		
 
  \item{filename}{the name of the file used to export the plot. \emph{character}.}
  \item{...}{additional arguments to be passed to \code{\link{optionsMRIaggr}} for specifying the graphical parameters.}
}
\details{
ARGUMENTS: \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

Information about the \code{lty}, \code{lwd}, \code{mar} and \code{mgp} arguments can be found in \code{\link{par}}.

Argument(s) \code{...} must correspond to some of the following arguments : \code{height}, \code{numeric2logical}, \code{path}, \code{res}, \code{unit}, \code{width}, \code{window}.

ARGUMENTS: \cr
Require to have previously filled \code{object@table_lesion}. This can be done using the \code{\link{calcTableLesion}} method.
}
\value{
  None.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## matplot display of the lesion
plotTableLesion(MRIaggr.Pat1_red, num = 1:3, type = "matplot",
                 mask = c("MASK_DWI_t0", "MASK_T2_FLAIR_t2"))

## evolution display of the lesion
plotTableLesion(MRIaggr.Pat1_red, num = 1:3, type = "evolution",
                mask = c("MASK_DWI_t0", "MASK_T2_FLAIR_t2"))

}
\concept{plot.}
\keyword{methods}

