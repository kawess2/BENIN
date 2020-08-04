\name{boxplotMask}
\title{Boxplot spatial group characteristics}
\alias{boxplotMask}
\alias{boxplotMask,MRIaggr-method}
\description{
Display a boxplot of the contrast parameter values inside and outside a spatial group.
}
\usage{
\S4method{boxplotMask}{MRIaggr}(object, param, mask, num = NULL, rescale = TRUE,  
         ylim = NULL, col = c("white", "purple"), main = NULL, mgp = c(2, 0.5, 0),
		 x.legend = "topright", y.legend = NULL, cex.legend = 0.8, 
         filename = paste(object@identifier, "boxplotMask", sep = "_"), ...)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
 
  \item{param}{the contrast parameter(s) associated with the lesion mask(s). \emph{character vector}. REQUIRED.}
  \item{mask}{the binary contrast parameter(s) defining the spatial group(s). \emph{character vector}. REQUIRED.}
  \item{num}{the slices to consider. \emph{numeric vector} or \code{NULL}.}
  \item{rescale}{should the contrast parameters be scaled ? \emph{logical}.}
   
  \item{ylim}{the y limits of the plot. \emph{numeric vector of size 2} or \code{NULL} leading to automatic setting of the y limits.}
  \item{col}{the colors of the boxplots for observations inside and outside the mask(s). \emph{character vector of size 2}.}
  \item{main}{an overall title for the plot. \emph{character}.}
  \item{mgp}{the margin line for the axis title, axis labels and axis line. \emph{positive numeric vector of size 3}.}
  \item{x.legend}{the x coordinates of the legend. \emph{numeric} or \emph{character}.}
  \item{y.legend}{the y coordinates of the legend. \emph{numeric} or \emph{character}.}
  \item{cex.legend}{the expansion factor of the legend. \emph{positive numeric}.}		
 
  \item{filename}{the name of the file used to export the plot. \emph{character}.}
  \item{...}{additional arguments to be passed to \code{\link{optionsMRIaggr}} for specifying the graphical parameters.}
}
\value{
  None.
}
\details{
ARGUMENTS: \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

Information about the \code{mgp} argument can be found in \code{\link{par}}.

Information about the \code{x.legend}, \code{y.legend}, \code{cex.legend} arguments can be found in \code{\link{legend}} (\code{cex.legend} is the \code{cex} argument of \code{legend}).

Arguments \code{...} must correspond to some of the following arguments : \code{height}, \code{hemisphere}, \code{norm_mu}, \code{norm_sigma}, \code{numeric2logical}, \code{path}, \code{res}, \code{unit}, \code{width}, \code{window}.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## display
boxplotMask(MRIaggr.Pat1_red, param = c("DWI_t0","TTP_t0","MTT_t0"),mask = "MASK_T2_FLAIR_t2",
     numeric2logical = TRUE)
}
\concept{plot.}
\keyword{methods}

