\name{calcROCthreshold}
\title{ROC analysis}
\alias{calcROCthreshold}
\alias{calcROCthreshold,MRIaggr-method}
\description{
Perform a ROC analysis on a continuous variable for identifying a binary outcome.
}
\usage{
\S4method{calcROCthreshold}{MRIaggr}(object, param, mask, plot = "ROC_Youden", digit = 10,
         filename = paste(object@identifier, "calcROCthreshold", plot, sep = "_"),
         update.object = FALSE, overwrite = FALSE, ...)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{param}{the contrast parameter(s) that should be used to identify the observations inside the mask. \emph{character vector}. REQUIRED.}
  \item{mask}{the binary contrast parameter that will be used as the outcome in the ROC analysis. \emph{character vector}. REQUIRED.}
  \item{plot}{the type of the graphic to display? \emph{character} or \code{FALSE}. See the details section.}
  \item{digit}{the number of digits used to round the data. \emph{positive integer}. \code{NULL} indicates no rounding.}
  \item{filename}{the name of the file used to export the plot. \emph{character}.}
  
  \item{update.object}{should the resulting threshold analysis be stored in \code{object@ls_descStats} ? \emph{logical}.}
  \item{overwrite}{if a threshold analysis is already stored in \code{object@ls_descStats}, can it be overwritten ? \emph{logical}.}  
  \item{...}{additional arguments to be passed to \code{\link{optionsMRIaggr}} for specifying the graphical parameters.}
}
\details{

This function requires to have installed the \emph{ROCR} package to work.

ARGUMENTS: \cr
Arguments \code{...} must correspond to some of the following arguments : \code{height}, \code{numeric2logical}, \code{path}, \code{res}, \code{unit}, \code{verbose}, \code{width}, \code{window}.

Possible values for \code{plot} are:
\itemize{
\item \code{"ROC_Youden"} : display the ROC curve with the optimal threshold according to the Youden index.
\item \code{"ROC_prev"} : display the ROC curve with the optimal threshold according to the utility function.
\item \code{"boxplot_Youden"} : display a boxplot of the contrast parameter for each outcome with the optimal threshold according to the Youden index.
\item \code{"boxplot_prev"} : display a boxplot of the contrast parameter for each outcome with the optimal threshold according to the utility function.
\item \code{FALSE} : no graphic is displayed.
}

}
\value{
   An \emph{data.frame} containing for each mask the AUC and AUPRC value, the optimal threshold and the corresponding sensitivity and specificity for the Youden criteria and a utility function taking into account the prevalence.
}
\seealso{
\code{\link{selectDescStats}} to select the mask caracteristics.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")
  
## ROC analysis
res <- calcROCthreshold(MRIaggr.Pat1_red, param = c("DWI_t0","T2_FLAIR_t2"),
                        mask = c("MASK_DWI_t0","MASK_T2_FLAIR_t2"), numeric2logical = TRUE)
  
res <- calcROCthreshold(MRIaggr.Pat1_red, param = c("DWI_t0","T2_FLAIR_t2"),
                        mask = c("MASK_DWI_t0","MASK_T2_FLAIR_t2"), numeric2logical = TRUE,
                        plot = "boxplot_Youden",
						update.object = TRUE, overwrite = TRUE)
						
selectDescStats(MRIaggr.Pat1_red, "Mask_threshold")
}
\concept{calc.}
\keyword{methods}

