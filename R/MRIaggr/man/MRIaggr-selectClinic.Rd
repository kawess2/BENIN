\name{selectClinic}
\title{Extract clinical data}
\alias{selectClinic}
\alias{selectClinic,MRIaggr-method}
\description{
Extract the clinical data from a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{selectClinic}{MRIaggr}(object, param = NULL)
}
\arguments{
  \item{object}{an \code{object} of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{param}{the clinical parameters to extract. \emph{character vector} or \code{NULL} leading to extract all the clinical parameters.}
}
\value{
  A one line \emph{data.frame} containing the clinical data in columns.
}
\seealso{
  \code{\link{allocClinic<-}} to allocate values in the clinic slot. 
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## select all clinical data
res <- selectClinic(MRIaggr.Pat1_red)

## select only the gender
res <- selectClinic(MRIaggr.Pat1_red, param = "Sex")
}
\concept{select.}
\keyword{methods}

