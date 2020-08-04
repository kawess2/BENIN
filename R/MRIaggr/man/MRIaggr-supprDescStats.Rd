\name{supprDescStats}
\title{Remove an element of ls_descStats}
\alias{supprDescStats<-}
\alias{supprDescStats<-,MRIaggr-method}
\description{
Remove an element from the \code{ls_descStat} attribute of a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{supprDescStats}{MRIaggr}(object, 
         verbose = optionsMRIaggr("verbose")) <- value
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{value}{the name of the element(s) that should be removed. \emph{character vector}. REQUIRED.}
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
}
\value{
  None.
}
\seealso{ 
 \code{\link{allocDescStats<-}} to allocate an element to the ls_descStat attribute. \cr
 \code{\link{selectParameter}} to display the elements of the ls_descStat attribute.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## existing elements in @ls_descStats
selectParameter(MRIaggr.Pat1_red, type = "ls_descStats")

## delete one element in @ls_descStats
supprDescStats(MRIaggr.Pat1_red) <- "index_sauve"

## remaining elements in @ls_descStats
selectParameter(MRIaggr.Pat1_red, "ls_descStats")
}
\concept{alloc.}
\keyword{methods}

