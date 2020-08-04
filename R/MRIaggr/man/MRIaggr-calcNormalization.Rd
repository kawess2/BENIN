\name{calcNormalization}
\title{Compute normalization values}
\alias{calcNormalization}
\alias{calcNormalization,MRIaggr-method}
\description{
Compute the normalization values for each contrast parameter.
}
\usage{
\S4method{calcNormalization}{MRIaggr}(object, param, mu_type = "mean", sigma_type = "sd",
         rm.CSF = FALSE, rm.WM = FALSE, rm.GM = FALSE,
         verbose = optionsMRIaggr("verbose"), update.object = FALSE, overwrite = FALSE)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{param}{the contrast parameters to normalize. \emph{character vector}. REQUIRED.}
  \item{mu_type}{the type of centering. Can be \code{"mean"} or \code{"median"}. }
  \item{sigma_type}{the type of scaling. Can be \code{"sd"} or \code{"mad"}.}
  \item{rm.CSF}{should the cerebral spinal fluid observations be excluded ? \emph{logical}.}
  \item{rm.GM}{should the grey matter observations be excluded ? \emph{logical}.}
  \item{rm.WM}{should the white matter observations be excluded ? \emph{logical}.}
 
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
  \item{update.object}{should the resulting normalization values be stored in \code{object} ? \emph{logical}.}
  \item{overwrite}{if normalization values are already stored in \code{object@normalization}, can they be overwritten ? \emph{logical}.}  
}
\details{
FUNCTION: \cr
If any of the \code{rm.CSF}, \code{rm.WM} or \code{rm.GM} is set to true, then the values of the parameters remaining to \code{FALSE} (among \code{CSF}, \code{WM} and \code{GM}) are summed. Voxels with value under 0.5 are discarded. \cr
Note that  \code{rm.CSF}, \code{rm.GM} and \code{rm.WM} cannot be set simultaneously to \code{TRUE}.
}
\value{
  An \emph{list} containing the normalization values, one element for each type of normalization.
}
\seealso{
\code{\link{selectNormalization}} to select the normalization values. \cr
\code{\link{calcTissueType}} to compute a probabilistic classification of the brain observations in WM/GM/CSF.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## compute normalization values
res <- calcNormalization(MRIaggr.Pat1_red, param = c("DWI_t0","T2_FLAIR_t2"),
                         update.object = TRUE, overwrite = TRUE)

## display
par(mfrow = c(2,4), mar = rep(1.5,4), mgp = c(2,0.5,0))
multiplot(MRIaggr.Pat1_red, param = "T2_FLAIR_t2", num = 1:3,
             legend = TRUE, window = NULL, main = "raw - slice ")
multiplot(MRIaggr.Pat1_red, param = "T2_FLAIR_t2", num = 1:3,
			 norm_mu="contralateral", norm_sigma="contralateral",
             legend = TRUE, window = NULL, main = "normalized - slice ")

## extract normalization
selectNormalization(MRIaggr.Pat1_red, type = "global", mu = TRUE, sigma = FALSE)

}
\concept{calc.}
\keyword{methods}

