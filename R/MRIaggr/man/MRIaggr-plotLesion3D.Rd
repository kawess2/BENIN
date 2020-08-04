\name{plotLesion3D}
\title{3D plot of the lesion}
\alias{plotLesion3D}
\alias{plotLesion3D,MRIaggr-method}
\description{
Make a 3D plot of the lesion. (experimental version)
}
\usage{
\S4method{plotLesion3D}{MRIaggr}(object, mask, edge = FALSE, Neighborhood = "3D_N6",
         numeric2logical = FALSE, spatial_res = c(1,1,1), xlim = NULL,
         ylim = NULL, zlim = NULL, type.plot = "shapelist3d", px_max = 10000,
         radius = 1, type = "s", col = "red", col.edge = "black")
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{mask}{the binary contrast parameter indicating the lesion. \emph{character}. REQUIRED.}
  \item{edge}{should the edges of the lesion be ploted instead of the core ? \emph{logical}.}
  \item{Neighborhood}{the type of neighbourhood used to defined the edges. \emph{character}.}
  \item{numeric2logical}{should \code{mask} be convert to logical ? \emph{logical}.}
  \item{spatial_res}{a dilatation factor for the coordinates. \emph{positive numeric vector of size 3}.}
  \item{xlim}{the x limits of the plot. \emph{numeric vector of size 2} or \code{NULL} leading to automatic setting of the x limits.}
  \item{ylim}{the y limits of the plot. \emph{numeric vector of size 2} or \code{NULL} leading to automatic setting of the y limits.}
  \item{zlim}{the z limits of the plot. \emph{numeric vector of size 2} or \code{NULL} leading to automatic setting of the z limits.}
  \item{type.plot}{the type of plot to be displayed. Can be \code{"plot3d"} or \code{"shapelist3d"}.}
  \item{px_max}{the maximum number of points that can be ploted. \emph{integer}.}
  \item{radius}{the radius of spheres. \emph{numeric}. See \code{plot3d} for more details.}
  \item{type}{the type of item to plot. \emph{character}. See \code{plot3d} for more details.}
  \item{col}{the color of the core of the lesion. \emph{character}.}
  \item{col.edge}{the color of the edge of the lesion. \emph{character}.}
}
\details{
ARGUMENTS: \cr
the \code{Neighborhood} argument can be a \emph{matrix} or an \emph{array} defining directly the neighbourhood to use (i.e the weight of each neighbor) 
or a name indicating which type of neighbourhood should be used (see the details section of \code{\link{initNeighborhood}}).

FUNCTION: \cr
This functions uses the \code{plot3d} or \code{shapelist3d} and thus require the \emph{rgl} package to work. This package is not automatically loaded by the \code{MRIaggr} package.
}
\value{
 None.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

\dontrun{
if(require(rgl)){
  # global view
plotLesion3D(MRIaggr.Pat1_red, mask = "MASK_T2_FLAIR_t2", spatial_res = c(1.875,1.875,6),
     numeric2logical = TRUE)

  # by slice
plotLesion3D(MRIaggr.Pat1_red, mask = "MASK_T2_FLAIR_t2", spatial_res = c(1.875,1.875,6),
     type.plot = "plot3d",
     numeric2logical = TRUE)
}
}
}
\concept{plot.}
\keyword{methods}

