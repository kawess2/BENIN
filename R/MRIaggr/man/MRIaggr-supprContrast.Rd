\name{supprContrast}
\title{Remove a contrast parameter}
\alias{supprContrast<-}
\alias{supprContrast<-,MRIaggr-method}
\description{
Remove a contrast parameter from a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{supprContrast}{MRIaggr}(object,
         verbose = optionsMRIaggr("verbose")) <- value
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{value}{the name of the parameter(s) that should be removed. \emph{character vector}. REQUIRED.}
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
}
\details{
ARGUMENTS: \cr
\code{value} can be a \emph{numeric vector} indicating the position of the parameters to remove in the \code{data} slot.
}
\seealso{ 
 \code{\link{allocContrast<-}} to allocate a contrast parameter. \cr
 \code{\link{selectParameter}} to display the contrast parameters.
}
\value{
  None.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## available contrast parameters
selectParameter(MRIaggr.Pat1_red)

## delete two contrast parameters
supprContrast(MRIaggr.Pat1_red) <- c("MTT_t0", "MASK_DWI_t0")

## remaining contrast parameters
selectParameter(MRIaggr.Pat1_red)
}
\concept{alloc.}
\keyword{methods}

