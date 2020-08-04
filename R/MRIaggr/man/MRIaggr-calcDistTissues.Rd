\name{calcDistTissues}
\title{Compute descriptive statistics}
\alias{calcDistTissues}
\alias{calcDistTissues,MRIaggr-method}
\description{
Compute the four first order statistics of the contrast parameters by cerebral structure.
}
\usage{
\S4method{calcDistTissues}{MRIaggr}(object, param, param.membership, num = NULL, hemisphere = "both")
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{param}{the contrast parameters to consider. \emph{character vector}. REQUIRED.}
  \item{param.membership}{the parameters indicating the probabilistic membership of each observations to each cerebral structure. \emph{character vector}.} 
  \item{num}{the slices to use. \emph{numeric vector} or \code{NULL}.}
  \item{hemisphere}{the hemisphere to use. \emph{character}.}
}
\details{
ARGUMENTS : \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

Information about the \code{hemisphere} argument can be found in the details section of \code{\link{selectContrast}}.
}
\value{
  An \emph{data.frame} containing in row the various cerebral structures and in columns the various moments for each contrast parameter.
}
\seealso{
\code{\link{calcTissueType}} to compute a probabilistic classification of the brain observations in WM/GM/CSF.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")
	   
## compute the distribution of DWI and T2 FLAIR for the CSF, WM, GM and lesion observations
res <- calcDistTissues(MRIaggr.Pat1_red, param = c("DWI_t0","T2_FLAIR_t2"),
                        param.membership = c("CSF","WM","GM","MASK_DWI_t0")
)
}
\concept{calc.}
\keyword{methods}

