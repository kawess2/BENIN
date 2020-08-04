\name{calcDistMask}
\title{Euclidean distance to a spatial group}
\alias{calcDistMask}
\alias{calcDistMask,MRIaggr-method}
\description{
Compute the euclidean distance to a spatial group.
}
\usage{
\S4method{calcDistMask}{MRIaggr}(object, mask, name_newparam = paste("dist", mask, sep = "_"),
     spatial_res = c(1,1,1), numeric2logical = FALSE, Neighborhood = "3D_N10",
	 verbose = optionsMRIaggr("verbose"), update.object = FALSE, overwrite = FALSE)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{mask}{the binary contrast parameter(s) defining the spatial groups from which the distance will be computed. \emph{character vector}. REQUIRED.}
  \item{name_newparam}{the name of the new distance parameters. \emph{character vector}.}
  \item{spatial_res}{a dilatation factor for the coordinates. \emph{positive numeric vector of size 3}.}
  \item{numeric2logical}{should \code{mask} be convert to logical ? \emph{logical}.}
  \item{Neighborhood}{the type of neighbourhood. \emph{character}.}
  
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.} 
  \item{update.object}{should the resulting distance parameters be stored in \code{object} ? \emph{logical}.}
  \item{overwrite}{if contrast parameters with the same names are already stored in \code{object@data}, can they be overwritten ? \emph{logical}.}
}  
\details{

This function requires to have installed the \emph{RANN} package to work.

ARGUMENTS: \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

The \code{Neighborhood} argument can be a \emph{matrix} or an \emph{array} defining directly the neighbourhood to use (i.e the weight of each neighbor) 
or a name indicating which type of neighbourhood should be used (see the details section of \code{\link{initNeighborhood}}).

FUNCTION: \cr
This function relies on the \code{nn2} function of the \emph{RANN} package.
}
\value{
  An \emph{data.frame} containing in row the observations and in columns the distance parameters.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## compute distance to initial and final mask
res <- calcDistMask(MRIaggr.Pat1_red, mask = c("MASK_DWI_t0","MASK_T2_FLAIR_t2"),
             update.object = TRUE, overwrite = TRUE)

multiplot(MRIaggr.Pat1_red, param = "dist_MASK_DWI_t0",
          index1 = list(coords = "MASK_DWI_t0", outline = TRUE))

## compute distance to initial and final mask correcting anisotropy
res <- calcDistMask(MRIaggr.Pat1_red, mask = c("MASK_DWI_t0","MASK_T2_FLAIR_t2"),
                    spatial_res = c(1.875,1.875,6),
                    update.object = TRUE, overwrite = TRUE)

multiplot(MRIaggr.Pat1_red, param = "dist_MASK_DWI_t0",
          index1 = list(coords = "MASK_DWI_t0", outline = TRUE))
}
\concept{calc.}
\keyword{methods}

