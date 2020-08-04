\name{selectDescStats}
\title{Extract non-stardard elements}
\alias{selectDescStats}
\alias{selectDescStats,MRIaggr-method}
\description{
Extract elements in the \code{ls_descStats} slot of a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{selectDescStats}{MRIaggr}(object, name = NULL)
}
\arguments{
  \item{object}{an \code{object} of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{name}{the name of the element to select. \emph{character} or \code{NULL} leading to select all available elements.}
}
\details{
This function requires to have installed the \emph{Matrix} and the \emph{spam} package to work when argument \code{name} is set to \code{"W_euclidean"}.
}
\seealso{ 
 \code{\link{allocDescStats<-}} to allocate elements in the \code{ls_descStats} slot. \cr
 \code{\link{calcHemisphere}} to identify the hemispheres. \cr
 \code{\link{allocHemisphere<-}} and \code{\link{allocContrast<-}} to allocate hemispheres obtained from an external source.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

calcGroupsMask(MRIaggr.Pat1_red,mask = c("MASK_DWI_t0","MASK_T2_FLAIR_t2"),
               W.range = 6, W.spatial_res = c(1.875,1.875,6),
               update.object = TRUE, overwrite = TRUE)

## select all elements in the slot @ls_descStats
ls_descStats <- selectDescStats(MRIaggr.Pat1_red)
names(ls_descStats)

## get the name of all elements present in the slot @ls_descStats
selectParameter(MRIaggr.Pat1_red, type = "ls_descStats")

## select a specific element 
res <- selectDescStats(MRIaggr.Pat1_red, name = "GroupsLesion")
}
\concept{select.}
\keyword{methods}

