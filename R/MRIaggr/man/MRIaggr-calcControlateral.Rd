\name{calcContralateral}
\title{Compute contraleral normalization values}
\alias{calcContralateral}
\alias{calcContralateral,MRIaggr-method}
\description{
Associate each voxel to an hemisphere and compute the difference between the voxel values and their contralateral correspondent.
}
\usage{
\S4method{calcContralateral}{MRIaggr}(object, param, num = NULL, type = "mean", 
     param.ref = NULL, distband = 1, lambda = 1, 
	 verbose = optionsMRIaggr("verbose"), update.object = FALSE, overwrite = FALSE)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{param}{the contrast parameters to normalize. \emph{character vector}. REQUIRED.}
  \item{num}{the slices to extract. \emph{numeric vector} or \code{NULL}.}
  \item{type}{the method used to compute the contralateral correspondent of each voxel. Can be \code{"mean"}, \code{"median"} or \code{"1NN_penalised"}.}
  \item{param.ref}{the parameter to use as a reference for the identification of the contralateral voxel. \emph{character} or \code{NULL} if no reference parameter available.}
  \item{distband}{the distance within which the contralateral values are considered. \emph{postive numeric}.}
  \item{lambda}{the importance of the penalization. \emph{numeric}.}
  
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
  \item{update.object}{should the resulting contralateral parameters be stored in \code{object} ? \emph{logical}.}
  \item{overwrite}{if contrast parameters with the same names are already stored in \code{object@data}, can they be overwritten ? \emph{logical}.}
}
\details{
ARGUMENTS: \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

FUNCTION: \cr
To compute the contralateral correspondent of each voxel, the mean or median value of the contralateral observations present in the distband can be used. 
Otherwise, considering a reference parameter, the contralateral voxel that minimised the difference in contrast (denoted I) with the voxel of interest (denoted x) penalized by the distance is retained :
\deqn{voxel_{contro} = argmin_{y} \frac{|I(x) - I(y)|}{\sigma} + \lambda * dist(x,y) where \sigma=sd(I)}
The \code{param.ref} and \code{lambda} are active only if \code{type} is \code{"1NN_penalized"}. In this case \code{lambda} equal \code{0} means no penalization.
}
\value{
  An \emph{list} containing :
  \itemize{
  \item\code{[[data]]} : a \emph{data.frame} containing the coordinates and the parameters normalized by the contralateral values.
  \item\code{[[index_plot]]} : two lists containing the observations used to compute the contralateral values, one for each hemisphere.
  }
}
\seealso{
\code{\link{calcHemisphere}} to identify the hemispheres.
\code{\link{selectContrast}} to select the contralateral normalized parameters or the hemisphere parameter.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## associate each voxel to its contralateral correspondant 
## according T1 parameter and compute the normalized parameters
res <- calcContralateral(MRIaggr.Pat1_red, param = c("DWI_t0","T2_FLAIR_t2"), num = NULL,
         type = "mean", param.ref = "T1_t0", distband = 1, lambda = 1, verbose = TRUE)

## display
multiplot(res$data[,c("i","j","k"),drop = FALSE],
             contrast = res$data$DWI_t0_contro
)

multiplot(res$data[,c("i","j","k"), drop = FALSE],
             contrast = res$data$DWI_t0_contro,
             index1 = res$data[res$index_plot$index.plot_lesionR,c("i","j","k"), drop = FALSE],
             index2 = res$data[res$index_plot$index.plot_lesionL,c("i","j","k"), drop = FALSE]
) 
}
\concept{calc.}
\keyword{methods}

