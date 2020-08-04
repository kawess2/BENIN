\name{calcBrainMask}
\title{Brain-Background discrimination}
\alias{calcBrainMask}
\alias{calcBrainMask,MRIaggr-method}
\description{
Seek to form two groups of observations (brain observations and background observations) using a threshold approach or a k-means algorithm.
}
\usage{
\S4method{calcBrainMask}{MRIaggr}(object, param, type = "kmeans", th.breaks = 100, 
         th.smoothing = TRUE, th.select_optima = 1, th.upper = TRUE, plot = TRUE,
         kmeans.n_groups = 2:4, kmeans.Neighborhood = 3, skull.param = NULL, 
		 skull.n_groups = 3,
		 filename = paste("calcBrainMask", type, object@identifier, sep = "_"),
         update.object = FALSE, overwrite = FALSE, ...)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{param}{the contrast parameter(s) that should be used to identify the brain observations. \emph{character vector}. REQUIRED.}
  \item{type}{the method to use. Can be \code{"threshold"} or \code{"kmeans"}.}
  
  \item{th.breaks}{the number of thresholds to use. \emph{postive integer}.}
  \item{th.smoothing}{should the derivative be smoothed ? \emph{logical}.}
  \item{th.select_optima}{the rank of the optimum to retain. \emph{postive integer}.}
  \item{th.upper}{should the observations above the selected threshold be retained ? Else the observations bellow will the selected threshold be retained. \emph{logical}.}
  \item{plot}{should the results be plotted ? \emph{logical}.}

  \item{kmeans.n_groups}{the number of groups to use in the kmeans algorithm. \emph{postive integer vector}.}
  \item{kmeans.Neighborhood}{the range of the neighbourhood. \emph{postive integer}.}
  \item{skull.param}{the parameter used to identify the skull. \emph{character}.}
  \item{skull.n_groups}{the number of groups to use in the kmeans algorithm to obtain the skull.}

  \item{filename}{the name of the file used to export the plot. \emph{character}.}

  \item{update.object}{should the resulting mask be stored in \code{object} as a \code{"mask"} parameter ? \emph{logical}.}
  \item{overwrite}{if a mask is already stored in \code{object@data}, can it be overwritten ? \emph{logical}.}
  \item{...}{additional arguments to be passed to \code{\link{optionsMRIaggr}} for specifying the graphical parameters.}

}
\details{
ARGUMENTS: \cr

Argument(s) \code{...} must correspond to some of the following arguments : \code{height}, \code{numeric2logical}, \code{path}, \code{res}, \code{unit}, \code{verbose}, \code{width}, \code{window}.
Setting \code{skull.param} to \code{NULL} leads to skip the skull stripping step.

FUNCTION: \cr
The threshold approach searches the best break point of the function that maps thresholds to the number of observations. For this, it find the optima of the first derivative of this function (possibly smoothed). \cr
The arguments \code{th.breaks}, \code{th.smoothing}, \code{th.select_optima} and \code{th.upper} are only active if \code{type} is \code{"threshold"}. \cr
\code{th.smoothing} can be set to an \emph{integer} to specify the width of the smoothing kernel.

The k-means approach seeks the most spatially coherent partition of the observations, among the possible partitions defined by the \code{kmeans.n_groups} argument. 
The mean number of neighbors averaged over observations (spatial potential) is used as a metric of the spatial coherence of the partition. \cr
The arguments \code{kmeans.n_groups}, \code{kmeans.Neighborhood} and \code{upper} are only active if \code{type} is \code{"kmeans"}.

The skull step consists in identifying the skull with an additionnal parameter (T1 sequence appears well suited for this purpose), and remove the corresponding observations from the brain mask. 
It is the partition that gives the best spatial coherence  for the final brain mask that is retained, leaving the posibility of no skull stripping.

The plot is active only if \code{type} is \code{"threshold"}.
}
\value{
If \code{type} is \code{"threshold"}, a \emph{list} containing :
\itemize{
\item\code{[[analysis]]} :  A \emph{matrix} containing the parameter thresholds (column \code{"threshold"}), the number of observations inside the mask (column \code{"Nb"}), the first and its smoothed version (column \code{"dNb"} and \code{"dNb.filtered"}) and indicator of optima (column \code{"optima"}). \emph{matrix}.
\item\code{[[th_opt]]} :  A \emph{matrix} containing the number of observations inside the mask and its derivative (in lines) for each optimum (in columns).
\item\code{[[best_group]]} :  An indicator variable giving the observations that belong to the mask. \emph{logical vector}.
\item\code{[[mask_name]]} :  the mask name. \emph{character}.
}         

If \code{type} is \code{"kmeans"}, a \emph{list} containing :
\itemize{
\item\code{[[kmeans]]} :  the optimal kmeans partition (result of the \code{kmeans} function). \emph{list}.
\item\code{[[potential]]} :  the spatial potential for the various brain partitions. \emph{matrix}.
\item\code{[[best_V]]} :  the highest potential. \emph{numeric}.
\item\code{[[best_group]]} :  An indicator variable giving the observations that belong to the mask. \emph{logical vector}.
\item\code{[[mask_name]]} : the mask name. \emph{character}.
\item\code{[[potential_skull]]} :  the spatial potential for the various skull partitions (only if \code{skull.param} is not \code{NULL}). \emph{matrix}.
}
} 
\seealso{
\code{\link{selectContrast}} to select the mask parameter. \cr
\code{\link{calcSmoothMask}} to spatially regularized the obtained mask.
}
\examples{
## load NIFTI files and convert them to MRIaggr
path.Pat1 <- system.file(file.path("nifti"), package = "MRIaggr")
ls.array <- list(readMRI(file.path(path.Pat1, "T1_t0"), format = "nifti"),
                 readMRI(file.path(path.Pat1, "T2_GRE_t0"), format = "nifti"))
MRIaggr.Pat1 <- constMRIaggr(ls.array, identifier = "Pat1", param = c("T1_t0","T2_GRE_t0"))

#### 1- thresholding approach ####
res <- calcBrainMask(MRIaggr.Pat1, param = "T2_GRE_t0", type = "threshold",
                     th.select_optima = 2)

breaks <- res$analysis[,"threshold"]
res <- calcBrainMask(MRIaggr.Pat1, param = "T2_GRE_t0", type = "threshold",
                     th.breaks = breaks[breaks > 50], th.select_optima = 1,
                     overwrite = TRUE, update.object = TRUE)

## display
multiplot(MRIaggr.Pat1, param = "mask")

multiplot(MRIaggr.Pat1, param = "T2_GRE_t0", index1 = "mask")

## other parameter 
\dontrun{
res <- calcBrainMask(MRIaggr.Pat1, param = "T1_t0", type = "threshold",
                     th.breaks = 200)

res <- calcBrainMask(MRIaggr.Pat1, param = "T1_t0", type = "threshold",
                     th.breaks = seq(0, 400, length.out = 50), th.select_optima = 2,
                     overwrite = TRUE, update.object = TRUE)

multiplot(MRIaggr.Pat1, param = "mask")
}

#### 2- k-means approach ####
\dontrun{
res <- calcBrainMask(MRIaggr.Pat1, param = "T2_GRE_t0", type = "kmeans",
                     kmeans.n_groups = 2:4,
                     update.object = TRUE, overwrite = TRUE)
}
\dontshow{
res <- calcBrainMask(MRIaggr.Pat1, param = "T2_GRE_t0", type = "kmeans",
                     kmeans.n_groups = 2,
                     update.object = TRUE, overwrite = TRUE)
}					 
## display
multiplot(MRIaggr.Pat1, param = "T2_GRE_t0", index1 = "mask")
multiplot(MRIaggr.Pat1, param = "mask")
}
\concept{calc.}
\keyword{methods}

