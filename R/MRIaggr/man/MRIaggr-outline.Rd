\name{outline}
\title{Outline a region on a slice}
\alias{outline}
\usage{
outline(n=50,sequential = TRUE, min_dist = 1,
         col = c("blue","red","grey"), pch = 20, cex = c(0.75,1,0.75))
}
\description{
Tool for graphical definition of a spatial region on an image.
}
\arguments{
  \item{n}{maximum number of points to define the outline. \emph{integer}.}
  \item{sequential}{should the region edge be updated on the graphical device after each point ? \emph{logical}.}
  \item{min_dist}{if the distance between the new point and the initial point is inferior to \code{min_dist}, then the definition of the region ends. \emph{numeric}. Only active if \code{sequential} is \code{TRUE}.}
  \item{col}{the colors in which the user-defined edge points, the interpolated edge points and the interior points should be ploted. \emph{character vector of size 3}.}
  \item{pch}{the symbol with which the observations will be displayed. \emph{positive integer}.}
  \item{cex}{the expansion factor used to plot the edge points, the interpolated edge points and the interior points. \emph{positive numeric vector of size 3}.}
}

\details{
FUNCTION: \cr
This function uses the \code{\link{locator}} function to obtain the coordinates of the cursor. It enable a point by point definition of a region where a linear interpolation is used between user-defined points to define the edge of the region. \cr
In the non sequential mode, the definition of the points stop if the number of points exceed \code{n} or using Echap. In the sequential mode, the definition of the points stop if the number of points exceed \code{n} or if the new point is close enough to the initial point. 

After defining the edge, the region is filled.
 }
\value{
  A list of two elements  :
   \itemize{
  \item\code{[[edge]]} : a \emph{data.frame} containing the position (\code{"i"} \code{"j"}) of the edge of the region, the edge number (\code{"edge"}), whether its user-specified point or interpolated point (\code{"points"}) and whether it was removed for the filling precedure (\code{"valid"}).
  \item\code{[[surface]]} : a \emph{data.frame} containing the position (\code{"i"} \code{"j"}) of the points belonging to the region.  
  }
}
\examples{
\dontrun{ 
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")
num <- 3

## display 1
multiplot(MRIaggr.Pat1_red, param = "T2_FLAIR_t2",              
             num = num, legend = FALSE, window = FALSE)

## outline on display 1
res <- outline(sequential=TRUE,min_dist=3)

## display the results
multiplot(MRIaggr.Pat1_red, param = "T2_FLAIR_t2",              
             num = num, legend = FALSE,
             index1 = data.frame(k = num, res$edge[,c("i","j")]),
             index2 = data.frame(k = num, res$surface[,c("i","j")]),
             window = FALSE)

carto <- selectContrast(MRIaggr.Pat1_red, param = c("MASK_T2_FLAIR_t2","index"),
             num = num, coords = TRUE)
carto <- merge(carto, cbind(res$surface, outline = TRUE), all = TRUE)
carto[is.na(carto$outline),"outline"] <- FALSE
head(carto)

## display the results next to MASK_T2_FLAIR_t2
multiplot(MRIaggr.Pat1_red, param = "T2_FLAIR_t2",              
             num = num, legend = FALSE,
             index1 = carto[carto$MASK_T2_FLAIR_t2,c("i","j","k")],
             index2 = carto[carto$outline,c("i","j","k")],
             window = FALSE)
}
}
\concept{calc.}
\keyword{functions}
