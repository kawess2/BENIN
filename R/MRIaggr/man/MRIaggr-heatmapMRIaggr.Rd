\name{heatmapMRIaggr}
\title{Correlation between contrast parameters}
\alias{heatmapMRIaggr}
\alias{heatmapMRIaggr,MRIaggr-method}
\description{
Display a correlation map of the contrast parameters.
}
\usage{
\S4method{heatmapMRIaggr}{MRIaggr}(object, param, num = NULL,
         rescale = TRUE, method = "pearson", points.values = TRUE, type = "image", 
		 breaks = NULL, col = cm.colors(256), main = NULL,
         mgp = c(2,0.5,0), mar = c(4,4,1,6), las = 1, cex.axis = 1,
         filename = paste(object@identifier, "heatmapMRIaggr", sep = "_"), ...)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{param}{the contrast parameters used to compute the correlations. \emph{character vector}. REQUIRED.}
  \item{num}{the slices to use. \emph{numeric vector} or \code{NULL}.}
  \item{rescale}{should the contrast parameters be scaled ? \emph{logical}.}

  \item{method}{the correlation coefficient which is to be computed. Can be \code{"pearson"}, \code{"kendall"} or \code{"spearman"}.}
  \item{points.values}{should the correlation values be printed on the plot ? \emph{logical}.}
  \item{type}{the type of plot to display. Any of \code{"image"} or \code{"image.plot"} or \code{FALSE} meaning no plot.}
  \item{breaks}{the break points to use to generate the color intervals. \emph{numeric vector} or \code{NULL} leading to automatic breakpoints generation.}

  \item{col}{the colors with which the correlations will be displayed. \emph{character vector}.}
  \item{main}{an overall title for the plot. \emph{character}.}
  \item{mgp}{the margin line for the axis title, axis labels and axis line. \emph{positive numeric vector of size 3}.}
  \item{mar}{the number of margin lines to be specified on the four sides of the plot. \emph{positive numeric vector of size 4}.}
  \item{las}{the style of the axis labels. Any of \code{0}, \code{1}, \code{2} or \code{3}.}
  \item{cex.axis}{the magnification to be used for axis annotation relative to the current setting of cex. \emph{positive numeric}.}

  \item{filename}{the name of the file used to export the plot. \emph{character}.}
  \item{...}{additional arguments to be passed to \code{\link{optionsMRIaggr}} for specifying the graphical parameters.}
}
\details{

This function with argument type set to \code{image.plot} requires to have installed the \emph{field} package to work.

Arguments \code{...} must correspond to some of the following arguments : \code{height}, \code{hemisphere}, \code{path}, \code{res}, \code{unit}, \code{width}, \code{window}.

ARGUMENTS: \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

If \code{breaks} is not \code{NULL}, it must be of length \code{length(col)+1}.

Information about the \code{mar}, \code{las} and \code{mgp} arguments can be found in \code{\link{par}}.
}
\value{
  None.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## pearson
heatmapMRIaggr(MRIaggr.Pat1_red, param = c("MASK_T2_FLAIR_t2","DWI_t0","TTP_t0","MTT_t0"),
           las = 1, type = "image", cex = 0.75,
           breaks = seq(-1, 1, length.out = 51),
           col = cm.colors(50))  	
## spearman
heatmapMRIaggr(MRIaggr.Pat1_red, param = c("MASK_T2_FLAIR_t2","DWI_t0","TTP_t0","MTT_t0"),
           las = 1, type = "image", cex = 0.75, method = "spearman",
           breaks = seq(-1, 1, length.out = 51),
           col = cm.colors(50))  


## spearman with legend
  heatmapMRIaggr(MRIaggr.Pat1_red, param = c("MASK_T2_FLAIR_t2","DWI_t0","TTP_t0","MTT_t0"),
           las = 1, type = "image.plot", cex = 0.75, method = "spearman",
           breaks = seq(-1, 1, length.out = 51),
           col = cm.colors(50))  

}
\concept{plot.}
\keyword{methods}

