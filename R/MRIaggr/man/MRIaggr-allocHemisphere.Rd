\name{allocHemisphere}
\title{Allocate the position of the mid-saggital plan}
\alias{allocHemisphere<-}
\alias{allocHemisphere<-,MRIaggr-method}
\usage{
\S4method{allocHemisphere}{MRIaggr}(object, 
         overwrite = FALSE, verbose = optionsMRIaggr("verbose")) <- value
}
\description{
Allocate information about the position of the mid-saggital plan to a \code{\linkS4class{MRIaggr}} object.
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{value}{a \emph{list} of \emph{data.frame}. Names must be among \code{"midplane"}, \code{"hemispheres"} \code{"data"}. See the Details section. REQUIRED.}  
  \item{overwrite}{if the characteristics of a mid-saggital plan are already stored in \code{object}, can they be overwritten ? \emph{logical}.}  
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
}
\details{
ARGUMENTS: \cr
The \code{"midplane"} element indicates the position of the observations of the mid-saggital plan. 
It has to be a two columns \emph{data.frame} with the coordinates \code{("i","j")} in columns and the observations in rows. 

The \code{"data"} element must contains the position of each voxel with respect to the mid-sagittal plan (column \code{"i_hemisphere"} and \code{"j_hemisphere"}) the hemisphere (\code{"left"}, \code{"right"} or \code{"undefined"}) to which the voxel belongs (column \code{"hemisphere"}). 

The \code{"hemispheres"} element indicates in which hemisphere is the lesion (denoted by \code{"lesion"}). The others hemispheries are denoted by \code{"contralateral"}.
It has to be a one line two columns \emph{data.frame} with names \code{"left"} \code{"right"}.
}
\seealso{
 \code{\link{calcHemisphere}} to obtain the position of the mid-saggital plan and the position of the lesion in the hemispheres. \cr
 \code{\link{selectParameter}} to extract the previous elements.
}
 \value{
  None.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## automatic allocation
resHemi <- calcHemisphere(MRIaggr.Pat1_red, param = "T2_GRE_t0", gridSearch = FALSE, num = 1,
                          update.object = TRUE, overwrite = TRUE)

## manual allocation
allocHemisphere(MRIaggr.Pat1_red, overwrite = TRUE) <- list(midplane = resHemi$midplane,
                                                            data = resHemi$data)

## display
index1 <- data.frame(selectMidplane(MRIaggr.Pat1_red), 15)
names(index1) <- c("i","j","k")
multiplot(MRIaggr.Pat1_red, param = "T2_GRE_t0", num = 1, midplane = TRUE, window = FALSE,
             index1 = list(coords = index1, pch = 20, cex = 2, col = "purple")
             
)

selectMidplane(MRIaggr.Pat1_red)
selectHemispheres(MRIaggr.Pat1_red)
}
\concept{alloc.}
\keyword{methods}

