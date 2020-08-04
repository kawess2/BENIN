\name{allocDescStats}
\title{Allocate non standard elements}
\alias{allocDescStats<-}
\alias{allocDescStats<-,MRIaggr-method}
\usage{
\S4method{allocDescStats}{MRIaggr}(object, name, 
         overwrite = FALSE, verbose = optionsMRIaggr("verbose")) <- value
}
\description{
Allocate non standard elements to a \code{\linkS4class{MRIaggr}} object.
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{value}{any R object. REQUIRED.}
  \item{name}{the name of the element storing \code{value}. \emph{character}. REQUIRED.}
  \item{overwrite}{if an element with the same name is already stored in \code{object@ls_descStats}, can it be overwritten ? \emph{logical}.}  
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
}
\details{
FUNCTION: \cr
Contrary to all other \code{alloc.} methods that impose restrictions on the objects that can be allocated, this function enable to allocate freely a R element.
However these elements will have no interaction with the methods of this package. 
}
\seealso{
 \code{\link{selectDescStats}} to extract non standard elements. 
}
\value{
  None.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## allocate a vector
allocDescStats(MRIaggr.Pat1_red,name = "spatial_res") <- c(1.875,1.875,6)

## select the corresponding element
selectDescStats(MRIaggr.Pat1_red, "spatial_res")

## some calc methods automatically save results in the ls_descStats slot
# find spatial groups 
calcGroupsMask(MRIaggr.Pat1_red,mask = c("MASK_DWI_t0","MASK_T2_FLAIR_t2"),
               W.range = 6, W.spatial_res = selectDescStats(MRIaggr.Pat1_red,"spatial_res"),
               update.object = TRUE, overwrite = TRUE)

# extract spatial groups
selectDescStats(MRIaggr.Pat1_red, "GroupsLesion")

}
\concept{alloc.}
\keyword{methods}
