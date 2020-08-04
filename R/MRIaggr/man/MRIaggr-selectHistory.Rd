\name{selectHistory}
\title{Extract the call of the methods applied on the object}
\alias{selectHistory}
\alias{selectHistory,MRIaggr-method}
\description{
Extract the history of a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{selectHistory}{MRIaggr}(object)
}
\arguments{
  \item{object}{an \code{object} of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
}
\value{
  A \emph{list} of each \code{calc} or \code{const} method that was applied to the object. 
  For each method, it contains the call, the date of call and potential extra elements.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## selection
selectHistory(MRIaggr.Pat1_red)
}        
\concept{select.}                                    
\keyword{methods}

