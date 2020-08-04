\name{selectTable}
\title{Extract volumic information}
\alias{selectTable}
\alias{selectTable,MRIaggr-method}
\description{
Extract the volumic information (lesion, hyperfusion or reperfusion) from a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{selectTable}{MRIaggr}(object, type, size = FALSE)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{type}{the table to extract. Can be \code{"lesion"}, \code{"reperfusion"} or \code{"hypoperfusion"}. REQUIRED.}
  \item{size}{should the values in the table correspond to a number of voxels (\code{FALSE}) or a volume (\code{TRUE}).}
}
\value{
  A \emph{data.frame} containing volumetric measurements.
}
\seealso{ 
 \code{\link{calcTableHypoReperf}} to compute and allocate the hypoperfusion and reperfusion tables. \cr 
 \code{\link{calcTableLesion}} to compute and allocate the lesion table. \cr
 \code{\link{allocTable<-}} to allocate data in the table slots. \cr 
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## selection
res <- selectTable(MRIaggr.Pat1_red, type = "lesion")
res <- selectTable(MRIaggr.Pat1_red, type = "reperfusion")
res <- selectTable(MRIaggr.Pat1_red, type = "hypoperfusion")
}
\concept{select.}
\keyword{methods}

