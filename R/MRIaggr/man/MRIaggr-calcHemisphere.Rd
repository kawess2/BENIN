\name{calcHemisphere}
\title{Find the mid-sagittal plan}
\alias{calcHemisphere}
\alias{calcHemisphere,MRIaggr-method}
\description{
Find a plane that distinguish the two cerebral hemispheres.
}
\usage{
\S4method{calcHemisphere}{MRIaggr}(object, param, num = NULL, p = 1, subset = NULL,
         penalty = "symmetry", mask = NULL, numeric2logical = FALSE,
		 n.points = 100, gridSearch = TRUE, 
         i_test = seq(-20, 20, by = 5), angle.test = seq(-30, 30, by = 5), 
		 unit_angle = "degree", NelderMead = TRUE, maxit = 100, reltol = 0.001, 
         plot = TRUE, filename = paste(object@identifier, "_calcHemisphere", sep = ""),
         update.object = FALSE, overwrite = FALSE, ...)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{param}{the contrast parameter that should be used to distinguish the two hemispheres. \emph{character}. REQUIRED.}
  \item{num}{the slices to use. \emph{numeric vector} or \code{NULL}.}
  \item{p}{the type of distance for the penalization. \emph{positive numeric}.}
  \item{subset}{the subset of observations to use. \emph{positive integer vector} or \code{NULL} leading to use all observations.}
  \item{penalty}{the type of objective function. Can be \code{"symmetry"} or \code{"asymmetry"}.}
  \item{mask}{the binary contrast parameter(s) indicating the lesion. \emph{character vector} or \code{NULL} if no mask is available.}
  \item{numeric2logical}{should \code{mask} be convert to logical ? \emph{logical}.}
 
  \item{gridSearch}{Should grid search be used to find the mid-sagittal plane ? \emph{logical}.}
  \item{i_test}{the abscissa or the number of abscissa to test. \emph{numeric vector} or \emph{positive integer}.}
  \item{angle.test}{the angle or the number of angle to test. \emph{numeric vector} or \emph{postive integer}.}
  \item{unit_angle}{the unit in which the angle is given. Can be \code{"radian"} or \code{"degree"}.}
  \item{n.points}{the number of points that represent the mid-sagittal plan to computed. \code{positive integer}.}

  \item{NelderMead}{Should the center of the grid search be searched using Nelder-Mead algorithm or set to the center of the image ? \emph{logical}.}
  \item{maxit}{The maximum number of iterations. \emph{postive integer}. See the details section of \code{optim}.}
  \item{reltol}{Relative convergence tolerance. \emph{positive numeric}. See the details section of \code{optim}.}
  
  \item{plot}{should the results be plotted ? \emph{logical}.}
  \item{filename}{the name of the file used to export the plot. \emph{character}.}
  \item{update.object}{should the resulting midplane be stored in \code{object} ? \emph{logical}.}
  \item{overwrite}{if a midplane is already stored in \code{object@midplane}, can it be overwritten ? \emph{logical}.}  
  \item{...}{additional arguments to be passed to \code{\link{optionsMRIaggr}} for specifying the graphical parameters.}
}
\details{
ARGUMENTS: \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

Setting \code{p} to \code{1} leads to use the absolute difference ; setting \code{p} to \code{2} leads to use the euclidean distance.

Arguments \code{...} must correspond to some of the following arguments : \code{height}, \code{numeric2logical}, \code{path}, \code{res}, \code{unit}, \code{verbose}, \code{width}, \code{window}.

FUNCTION: \cr
This function seeks the plane that minimize the difference between contralateral values of the two hemispheres or maximize the similarity between the two hemispheres.
There are 2 degree of freedom : one for the position of the center and one for the angle.
The separation between the hemisphere is assumed to be identical for all slices. 
From our experience, using an objective function based on symmetry gives better results compared to asymmetry.
}
\seealso{
\code{\link{selectParameter}} to select the midplane.
}
\value{
 An \emph{list} containing :
  \itemize{
  \item\code{[[penalty]]} : an \emph{array} containing the objective function for the various configurations.
  \item\code{[[nb]]} : an \emph{array} containing the number of observations used to compute the penalty function for the various configurations.
  \item\code{[[moy]]} : an \emph{array} containing the mean value of the objective function for the various configurations.
  \item\code{[[optimum]]} : the parameters of the optimal midplane. \emph{numeric vector}
  \item\code{[[midplane]]} : the position of the midplane points. \emph{matrix}.
  \item\code{[[data]]} : the position of the observations with respect to the mid-sagittal plane. \emph{matrix}.
 \item\code{[[cv]]} : Was the optimum reached inside the parameter space and not at a border ? \emph{logical}.
  }
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

\dontshow{
 res <- calcHemisphere(MRIaggr.Pat1_red, param = "T2_GRE_t0", num = 1,
                        mask = c("MASK_DWI_t0","MASK_T2_FLAIR_t2"), i_test = c(-5,0,5), angle.test = c(-0.5,0,0.5),
                        verbose = TRUE, update.object = TRUE, overwrite = TRUE)
}

\dontrun{
res <- calcHemisphere(MRIaggr.Pat1_red, param = "T2_GRE_t0",
                      verbose = TRUE, update.object = TRUE, overwrite = TRUE)

## display the mid-sagittal plan
multiplot(MRIaggr.Pat1_red, param = "T2_GRE_t0", num = 3, legend = FALSE,
             midplane = TRUE, main = "original coordinates - slice ")

## display with the new coordinates
multiplot(selectContrast(MRIaggr.Pat1_red, param = c("i_hemisphere","j_hemisphere","k")),
          contrast=selectContrast(MRIaggr.Pat1_red, param = "T2_GRE_t0"), num = 3,
          index1=cbind(0, seq(-50,50), 3), main = "new coordinates - slice ", legend = FALSE)

## compute the mid-sagittal plan and mark lesion and healthy hemispheres
res <- calcHemisphere(MRIaggr.Pat1_red, param = "T2_GRE_t0",
         mask = c("MASK_DWI_t0","MASK_T2_FLAIR_t2"), numeric2logical = TRUE,
         verbose = TRUE, update.object = TRUE, overwrite = TRUE)
}

}
\concept{calc.}
\keyword{methods}

