\name{calcSmoothMask}
\title{Spatial regularization}
\alias{calcSmoothMask}
\alias{calcSmoothMask,MRIaggr-method}
\description{
Perform a spatial regularization of a binary mask.
}
\usage{
\S4method{calcSmoothMask}{MRIaggr}(object, mask = "mask", numeric2logical = FALSE,
         size_2Dgroup = 50, Neighborhood_2D = "3D_N8", rm.2Dhole = FALSE,
         size_3Dgroup = "unique", Neighborhood_3D = "3D_N10", rm.3Dhole = TRUE,
         erosion.th = 0.75, Vmask_min = 0.25, Vbackground_max = 0.75,
         Neighborhood_V = "3D_N10", verbose = optionsMRIaggr("verbose"), 
		 update.object = FALSE, overwrite = FALSE)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{mask}{the binary contrast parameter that should be smoothed. \emph{character}.}
  \item{numeric2logical}{should \code{mask} be convert to logical ? \emph{logical}.}
 
  \item{size_2Dgroup}{the minimum size of the 2D groups. \emph{positive integer} or \code{"unique"}.}
  \item{Neighborhood_2D}{the type of 2D neighbourhood. \emph{character}.}
  \item{rm.2Dhole}{should the 2D wholes inside the mask be removed ? \emph{logical}.}
 
  \item{size_3Dgroup}{the minimum size of the 3D groups. \emph{positive integer} or \code{"unique"}.}
  \item{Neighborhood_3D}{the type of 3D neighbourhood. \emph{character}.}
  \item{rm.3Dhole}{should the 3D wholes inside the mask be removed ? \emph{logical}.}
  \item{erosion.th}{the threshold below which the observations will be removed by the erosion. \emph{numeric between 0 and 1}.}
 
  \item{Vmask_min}{mask observations with a proportion of neighbors belonging to the mask lower than \code{Vmask_min} are attributed to the background. \emph{numeric between 0 and 1}.}
  \item{Vbackground_max}{background observations with a proportion of neighbors belonging to the mask higher than \code{Vbackground_max} are attributed to the mask. \emph{numeric between 0 and 1}.}
  \item{Neighborhood_V}{the type of neighbourhood to use for the spatial regularization. \emph{character}.}
  
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
  \item{update.object}{should the resulting regularized mask be stored in \code{object} ? \emph{logical}.}
  \item{overwrite}{if a mask is already stored in \code{object@data}, can it be overwritten ? \emph{logical}.}  
}
\details{
ARGUMENTS: \cr
the \code{Neighborhood_2D} or \code{Neighborhood_3D} arguments can be a \emph{matrix} or an \emph{array} defining directly the neighbourhood to use (i.e the weight of each neighbor) 
or a name indicating which type of neighbourhood should be used (see the details section of \code{\link{initNeighborhood}}).

FUNCTION: \cr
This function applies 6 smoothing steps : 
\itemize{
\item exclusion of the small 2D groups from the mask (to skip set \code{size_2Dgroup} to \code{FALSE}). Note that \code{size_2Dgroup = "unique"} lead to keep the largest 2D group of each slice.
\item filling of the small 2D holes in the mask (to skip set \code{rm.2Dhole} to \code{FALSE}).
\item exclusion of the small 3D groups from the mask (to skip set \code{size_3Dgroup} to \code{FALSE}). Note that \code{size_3Dgroup = "unique"} lead to keep only the largest 3D group.
\item erosion that first temporarily remove observations from the mask that have less than \code{erosion.th} percent of their neighbourhood in the mask. 
Then it computes the new 3D groups and remove permanently all the new 3D groups from the mask. To skip set \code{erosion.th} to \code{FALSE}.
\item filling of the small 3D holes in the mask (to skip set \code{rm.3Dhole} to \code{FALSE}).
\item spatial regularization that homogenize the local neighbourhood (to skip set both \code{Vmask_min} and \code{Vbackground_max} to \code{FALSE}). 
}
}
\seealso{
\code{\link{selectContrast}} to select the smoothed mask.
\code{\link{calcBrainMask}} to compute an indicator of the brain observations.
}
\value{
  An \emph{data.frame} containing the mask and the coordinates in columns.
}
\examples{
## load data and build MRIaggr
path.Pat1 <- system.file(file.path("nifti"), package = "MRIaggr")
ls.array <- list(readMRI(file.path(path.Pat1,"T2_GRE_t0"), format = "nifti"))
MRIaggr.Pat1 <- constMRIaggr(ls.array,identifier="Pat1", param = "T2_GRE_t0")

## create the cerebral mask
res <- calcBrainMask(MRIaggr.Pat1, param = "T2_GRE_t0", type = "kmeans",
                     kmeans.n_groups = 2:4,
                     update.object = TRUE, overwrite = TRUE)

## smooth the cerebral mask
res <- calcSmoothMask(MRIaggr.Pat1, update.object = TRUE, overwrite = TRUE)

## display
multiplot(MRIaggr.Pat1,param = "mask", legend = FALSE)
}
\concept{calc.}
\keyword{methods}

