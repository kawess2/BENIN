\name{selectNormalization}
\title{Extract the normalization values}
\alias{selectNormalization}
\alias{selectNormalization,MRIaggr-method}
\description{
Extract the normalization values from a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{selectNormalization}{MRIaggr}(object, type = NULL, mu = TRUE, sigma = TRUE,
         hemisphere = "both", num = NULL, param = NULL)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{type}{the type of normalization. Must be on of \code{"global"}, \code{"slice"}, \code{"3slices"} or \code{NULL} leading to select all types of normalization.}
  \item{mu}{should the centering values for the normalization be returned. \emph{logical}.  Active only if \code{type} is \code{"slice"} or \code{"3slices"}.}
  \item{sigma}{should the scaling values for the normalization be returned. \emph{logical}.  Active only if \code{type} is \code{"slice"} or \code{"3slices"}.}
  \item{num}{the slices to extract. \emph{numeric vector} or \code{NULL}.}
  \item{hemisphere}{the hemisphere to extract. \emph{character}.}
  \item{param}{the contrast parameters for which the normalization values should be extracted. \emph{character vector} or \code{NULL} indicating all available contrast parameters.}
}
\details{
ARGUMENTS: \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

Information about the \code{hemisphere} argument can be found in the details section of \code{\link{selectContrast}}.
}
\value{
  A \emph{data.frame} or a \emph{list}.
}
\seealso{
 \code{\link{calcNormalization}} to compute and allocate the normalisation values. \cr
 \code{\link{allocNormalization<-}} to alloacte the normalization values when obtained from an external source. \cr
 \code{\link{calcHemisphere}} to identify the hemispheres. \cr
 \code{\link{allocHemisphere<-}} and \code{\link{allocContrast<-}} to allocate hemispheres obtained from an external source.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## select all normalization values
res <- selectNormalization(MRIaggr.Pat1_red)
names(res)

## select specific normalization normalization values 
# computed on the whole brain
res <- selectNormalization(MRIaggr.Pat1_red, type = "global", 
         mu = TRUE, sigma = FALSE, hemisphere = "both")
# idem but only for DWI_t0
res <- selectNormalization(MRIaggr.Pat1_red, type = "global",
         mu = TRUE, sigma = FALSE, param = "DWI_t0") 

# computed by slice
res <- selectNormalization(MRIaggr.Pat1_red, type="slice", 
         mu = TRUE,sigma = FALSE, hemisphere = "both") 
# idem for slice 1
res <- selectNormalization(MRIaggr.Pat1_red, type = "slice", 
         mu = TRUE, sigma = FALSE, num = 1) 

# computed on 3 consecutive slices
res <- selectNormalization(MRIaggr.Pat1_red, type = "3slices", mu = FALSE, sigma = TRUE, 
                           hemisphere = "left", num = 2, param = "T2_FLAIR_t2")
}
\concept{select.}
\keyword{methods}

