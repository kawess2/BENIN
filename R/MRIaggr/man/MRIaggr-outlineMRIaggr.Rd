\name{outlineMRIaggr}
\title{Outline a region on a slice}
\alias{outlineMRIaggr}
\alias{outlineMRIaggr,MRIaggr-method}
\usage{
\S4method{outlineMRIaggr}{MRIaggr}(object, param, index1 = NULL, num = NULL, hemisphere = "both",
     numeric2logical = FALSE, xlim = NULL, ylim = NULL, legend = FALSE,
     palette = "terrain.colors", col = NULL, breaks = 25, fill = TRUE, n = 50,
     sequential = TRUE, min_dist = 1, operator_index1 = "none",
     col.outline = c("blue","red","grey"), pch = 20, cex = c(0.75,1,0.75), 
	 name_newparam = "userMask",
	 verbose = optionsMRIaggr("verbose"), update.object = FALSE, overwrite = FALSE)
}
\description{
Tool for graphical definition of a spatial region on an image.
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{param}{the contrast parameter map on which the outline will be drawn. \emph{character}. REQUIRED.}
  \item{index1}{the coordinates of additionnal points to display. \emph{data.frame} or \emph{list} or \code{NULL}.}
  \item{num}{the slices on which the outline will be drawn. \emph{numeric vector} or \code{NULL}.}
  \item{hemisphere}{the hemisphere to consider. \emph{character}.}
  \item{numeric2logical}{if a parameter is specified for the index argument, should it be converted to logical ? \emph{logical}.}
 
  \item{xlim}{the y limits of the plot. \emph{numeric vector of size 2} or \code{NULL} leading to automatic setting of the x limits.}
  \item{ylim}{the y limits of the plot. \emph{numeric vector of size 2} or \code{NULL} leading to automatic setting of the y limits.}
  \item{legend}{how should the legend be displayed ? \emph{logical} or \code{NULL}.}
  \item{palette}{the colors or the palette to use when associating colors to observation intensities. \emph{character vector} or \emph{character}.}
  \item{col}{the color to use to plot the observations. \emph{character vector} or \code{NULL} leading to automatic generation of the colors using the \code{breaks} and \code{palette} arguments.}
  \item{breaks}{the break points or the number of breakpoints to use to generate the color intervals . \emph{numeric vector} or \emph{postive integer}.}
 
  \item{fill}{should the spatial region be filled ? Otherwise only the edge is used. \emph{logical}.}
  \item{n}{maximum number of points to define the outline. \emph{integer}.}
  \item{sequential}{should the region edge be updated on the graphical device after each point ? \emph{logical}.}
  \item{min_dist}{if the distance between the new point and the initial point is inferior to \code{min_dist}, then the definition of the region ends. \emph{numeric}. Only active if \code{sequential} is \code{TRUE}.}
  \item{operator_index1}{the operator to apply between the index1 observations and the outlined observations. Can be \code{"none"} \code{"intersection"} \code{"difference"} or \code{"union"}.}

  \item{col.outline}{the colors in which the user-defined edge points, the interpolated edge points and the interior points should be ploted. \emph{character vector}[3].}
  \item{pch}{the symbol with which the observations will be displayed. \emph{positive integer}.}
  \item{cex}{the expansion factor used to plot the edge points, the interpolated edge points and the interior points. \emph{positive numeric vector of size 3}.}
 
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
  \item{name_newparam}{the name of the new parameter. \emph{character}.}
  \item{update.object}{should the new parameter be stored in \code{object} ? \emph{logical}.}
  \item{overwrite}{if a contrast parameter with the same names is already stored in \code{object@data}, can it be overwritten ? \emph{logical}.}  
 }

\details{
ARGUMENTS: \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

Information about the \code{index1} argument can be found in the details section of \code{\link{initIndex}} (argument \code{index}).

Information about the \code{hemisphere} argument can be found in the details section of \code{\link{selectContrast}}.

Information about the \code{palette} and \code{legend} arguments can be found in the details section of \code{\link{initCol}}.

FUNCTION: \cr
This function uses the \code{\link{locator}} function to obtain the coordinates of the cursor. See \code{\link{outline}} for more details. 
It will display slice by slice the parameter map and on each slice a region should be drawn. \cr
Press Echap and 1 to restart the draw on the same slice. \cr
Press Echap and 0 to skip a slice. 
}
\value{
  A matrix containing in columns :
 \itemize{  
  \item\code{[,c("i","j","k")]} : the coordinates of the observations.
  \item\code{[,"index"]} : the index of the observation in \code{object}.
  \item\code{[,"edge"]} : whether each observation is on an edge.
  \item\code{[,"surface"]} : whether each observation is in the interior.
  \item\code{[,"userMask"]} : whether each observation belong to \code{index1} (if no not defined, always \code{FALSE}).
  } 
}
\examples{
\dontrun{ 
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## outline the area of interest
res <- outlineMRIaggr(MRIaggr.Pat1_red, param = "DWI_t0",
             num = 1:3, sequential = TRUE, overwrite = TRUE, update.object = TRUE)

multiplot(MRIaggr.Pat1_red, param = "userMask",              
             num = 1:3)


## outline an edge of interest 
res <- outlineMRIaggr(MRIaggr.Pat1_red, param = "DWI_t0", index1 = "MASK_DWI_t0",
            fill = FALSE, num = 1:3, sequential = TRUE, overwrite = TRUE, update.object = TRUE)

multiplot(MRIaggr.Pat1_red, param = "userMask", num = 1:3)

## define an new area as the union of the outlined area and the initial lesion mask
res <- outlineMRIaggr(MRIaggr.Pat1_red, param = "DWI_t0", index1 = "MASK_DWI_t0",
            operator_index1 = "union", num = 3, sequential = TRUE, 
			overwrite = TRUE, update.object = TRUE)

multiplot(MRIaggr.Pat1_red, param = "userMask", num = 3)

## define an new area as the intersection of the outlined area and the initial lesion mask
res <- outlineMRIaggr(MRIaggr.Pat1_red, param = "DWI_t0", index1 = "MASK_DWI_t0",
            operator_index1 = "intersection", num = 3, sequential = TRUE, 
			overwrite = TRUE, update.object = TRUE)

multiplot(MRIaggr.Pat1_red, param = "userMask", num = 3)

## define an new area as the difference between the outlined area and the initial lesion mask
res <- outlineMRIaggr(MRIaggr.Pat1_red, param = "DWI_t0", index1 = "MASK_DWI_t0",
            operator_index1 = "difference", num = 3, sequential = TRUE,
			overwrite = TRUE, update.object = TRUE)

multiplot(MRIaggr.Pat1_red, param = "userMask", num=3)
}
}
\concept{calc.}
\keyword{methods}
