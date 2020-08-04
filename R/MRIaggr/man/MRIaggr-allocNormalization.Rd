\name{allocNormalization}
\title{Allocate normalization values}
\alias{allocNormalization<-}
\alias{allocNormalization<-,MRIaggr-method}
\usage{
\S4method{allocNormalization}{MRIaggr}(object, 
         overwrite = FALSE, verbose = optionsMRIaggr("verbose")) <- value
}
\description{
Allocate normalization values to a \code{\linkS4class{MRIaggr}} object.
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{value}{the normalisation values. A \emph{list} of \emph{data.frame}. REQUIRED.}
  \item{overwrite}{if normalization values are already stored in \code{object}, can they be overwritten ? \emph{logical}.}  
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
}
\details{
ARGUMENTS: \cr
To enable other methods of the package to use the normalization values the \code{value} argument should match the result of the \code{calcNormalization} function.
This function only partially control the validity of the \code{value} argument. \cr
A valid \code{value} argument should be composed a list of \emph{data.frame} with the following names :
\itemize{
\item \code{"norm_global"} : it should have 6 rows named \code{"mu_both"} \code{"mu_left"} \code{"mu_right"} \code{"sigma_both"} \code{"sigma_left"} \code{"sigma_right"}.
\item \code{"normMu_slice_both"} :
\item \code{"normSigma_slice_both"}
\item \code{"normMu_slice_left"} 
\item \code{"normSigma_slice_left"}   
\item \code{"normMu_slice_right"} 
\item \code{"normSigma_slice_right"} 
\item \code{"normMu_3slices_both"} 
\item \code{"normSigma_3slices_both"} 
\item \code{"normMu_3slices_left"} 
\item \code{"normSigma_3slices_left"} 
\item \code{"normMu_3slices_right"} 
\item \code{"normSigma_3slices_right"}.
}
All the elements of the list must have column names that match the contrast parameters present in the \code{object}. \cr 
Appart from the first element of the list, all elements should have as many rows as slices contained in the object. 
}
\seealso{
 \code{\link{calcNormalization}} to compute the normalisation values. \cr
 \code{\link{selectNormalization}} to extract the normalisation values. 
}
\value{
  None.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## parameters to normalize
param <- c("DWI_t0","T2_FLAIR_t2","T2_GRE_t0","TTP_t0")

## manual allocation
resNormalization <- calcNormalization(MRIaggr.Pat1_red, param = param)
allocNormalization(MRIaggr.Pat1_red, overwrite = TRUE) <- resNormalization

## automatic allocation
resNormalization <- calcNormalization(MRIaggr.Pat1_red, param = param,
                                      update.object = TRUE, overwrite = TRUE)

}
\concept{alloc.}
\keyword{methods}

