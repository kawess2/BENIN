\name{calcThreshold}
\title{Image thresholding}
\alias{calcThreshold}
\alias{calcThresholdMRIaggr}
\alias{calcThresholdMRIaggr,MRIaggr-method}
\description{
Threshold a contrast parameter at one or several values.
}
\usage{
\S4method{calcThresholdMRIaggr}{MRIaggr}(object, param, hemisphere = "both", rm.CSF = FALSE,
     threshold = 1:10, decreasing = FALSE,
	 GRalgo = FALSE, W = "ifany", seed = NULL, numeric2logical = FALSE, W.range,
	 W.spatial_res = rep(1,3), name_newparam = paste(param,"Th", sep = "_"),
	 verbose = optionsMRIaggr("verbose"), update.object = FALSE, overwrite = FALSE)
	 
calcThreshold(contrast, param, hemisphere = NULL, rm.CSF = FALSE, threshold = 1:10,
     decreasing = FALSE, GRalgo = FALSE, W = NULL, seed = NULL,
     numeric2logical = FALSE, verbose = optionsMRIaggr("verbose"))
}				
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{contrast}{the dataset containing the contrast parameter to be thresholded. \emph{matrix}. REQUIRED.}
  \item{param}{the contrast parameters to be thresholded. \emph{character vector}. REQUIRED.}
  \item{hemisphere}{the hemisphere to consider. \emph{character} or \code{NULL}.}
  \item{rm.CSF}{should the cerebral spinal fluid observations be excluded ? \emph{logical} or \emph{character}.}
  \item{threshold}{the thresholds to be used for the discretization of the contrast parameter. \emph{numeric vector}.}
  \item{decreasing}{should the increasing thresholding (\code{FALSE}) or decreasing thresholding (\code{TRUE}) be used. \emph{logical}.}
  
  \item{GRalgo}{should a Growing Region algorithm be used to clean the thresholded parameter ? \emph{logical}.}
  \item{W}{the neighbourhood matrix. \emph{dgCMatrix} or \code{"ifany"} leading to use the neighbourhood matrix stored in the \code{object} if any and else computate this matrix.}
  \item{seed}{the index of the seeds for the growing region algorithm . \emph{positive integer vector}. }
  \item{numeric2logical}{should \code{seed} be converted to logical ? \emph{logical}.}
  \item{W.range}{only distances smaller than \code{W.range} are stored in \code{W}. \emph{postive numeric}.}
  \item{W.spatial_res}{a dilatation factor for the coordinates. \emph{positive numeric vector of size 3}.}
    
  \item{name_newparam}{the name of the new parameters. \emph{character}.}
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
  \item{update.object}{should the resulting thresholded parameters be stored in \code{object} ? \emph{logical}.}
  \item{overwrite}{contrast parameters with the same names are already stored in \code{object@data}, can it be overwritten ? \emph{logical}.}  
}
\details{
These functions require to have installed the \emph{Matrix} and the \emph{spam} package to work when argument \code{GRalgo} is set to \code{TRUE}.

ARGUMENTS: \cr
[data.frame method] the \code{hemisphere} argument must be one of the levels present in the column named \code{"hemisphere"} in data. \code{NULL} leading to use all observations. \cr
[MRIaggr method] Information about the \code{hemisphere} arguments can be found in the details section of \code{\link{selectContrast}}. 

FUNCTION: \cr
[data.frame method] By default the indicator of CSF will be extract from the column named \code{CSF}. 
If it is contained in another column the user must specify \code{rm.CSF} with the name of this column. \cr
[MRIaggr method] Setting \code{rm.CSF} to \code{TRUE} require to have the corresponding parameter (by default \code{CSF}) stored in the \code{object}.
It can be done using \code{calcTissueType}.
}
\seealso{
\code{\link{selectContrast}} to select contrast parameters in the \code{MRIaggr} object. 
}
\examples{
## load a MRIaggr object
data(MRIaggr.Pat1_red, package = "MRIaggr")

#### 1- MRIaggr method 
## raw parameter
multiplot(MRIaggr.Pat1_red, param = "TTP_t0", legend = FALSE, main = "TTP_t0 - slice ",
             palette=rainbow(10), breaks = seq(0,10) - 10^{-10})

## thresholded parameter
calcThresholdMRIaggr(MRIaggr.Pat1_red, param = c("TTP_t0","MTT_t0"),threshold = 1:10,
                 name_newparam = c("TTP.th_t0","MTT.th_t0"),
                 update.object = TRUE, overwrite = TRUE)

multiplot(MRIaggr.Pat1_red, param = "TTP.th_t0", main = "TTP.th_t0 - slice",
             legend = FALSE,palette = rainbow(10), breaks = (0:10) - 10^{-10})

\dontrun{
## 1st correction
calcThresholdMRIaggr(MRIaggr.Pat1_red, param = c("TTP_t0","MTT_t0"),threshold = 1:10,
                 rm.CSF = TRUE, hemisphere = "lesion", 
				 name_newparam = c("TTP.red_t0","MTT.red_t0"),
				 update.object = TRUE, overwrite = TRUE)

multiplot(MRIaggr.Pat1_red, param = "TTP.red_t0", main = "TTP.red_t0 - slice",
             legend = FALSE, palette = rainbow(10), breaks = (0:10) - 10^{-10})

## 2nd correction
calcThresholdMRIaggr(MRIaggr.Pat1_red, param = c("TTP_t0","MTT_t0"), threshold = 1:10,
         rm.CSF = TRUE, hemisphere = "lesion", name_newparam = c("TTP.GR_t0","MTT.GR_t0"),
         GRalgo = TRUE, seed = c("MASK_T2_FLAIR_t2","MASK_DWI_t0"), W.range  = sqrt(2),
         update.object = TRUE, overwrite = TRUE)

multiplot(MRIaggr.Pat1_red, param = "TTP.GR_t0", main = "TTP.GR_t0 - slice",
             legend = FALSE, palette = rainbow(10), breaks = (0:10) - 10^{-10})


#### 2- data.frame function ####

## raw parameter
multiplot(MRIaggr.Pat1_red, param = "TTP_t0", legend = FALSE, main = "TTP_t0 - slice ",
         palette=rainbow(10), breaks = seq(0,10) - 10^{-10})

## thresholded parameter 
data <- selectContrast(MRIaggr.Pat1_red, 
                       param = c("TTP_t0","MTT_t0","hemisphere","CSF","WM","GM"))
hypo_Th1_10 <- calcThreshold(data,param = c("TTP_t0","MTT_t0"), threshold = 1:10)

multiplot(selectCoords(MRIaggr.Pat1_red), main = "TTP_t0_th - slice ",
             hypo_Th1_10[,1], legend = FALSE, palette = rainbow(10), breaks = (0:10) - 10^{-10})

## 1st correction
data$CSF <- as.numeric(apply(data[,c("CSF","WM","GM")], 1, which.max) == 1)

hypoC_Th1_10 <- calcThreshold(data,param = c("TTP_t0","MTT_t0"), threshold = 1:10,
                              hemisphere = "left", rm.CSF = TRUE)

multiplot(selectCoords(MRIaggr.Pat1_red), main = "TTP_t0_thC - slice",
          hypoC_Th1_10[,1], legend = FALSE, palette = rainbow(10), breaks = (0:10) - 10^{-10})

## 2nd correction
maskN <- c("MASK_T2_FLAIR_t2","MASK_DWI_t0")
data[, maskN] <- selectContrast(MRIaggr.Pat1_red, param = maskN)
W <- calcW(MRIaggr.Pat1_red, range  = sqrt(2 * 1.875^2 + 0.001), row.norm = TRUE, upper = NULL,
           spatial_res = c(1.875,1.875,6))$W
max(spam::rowSums(W > 0))

hypoCC_Th1_10 <- calcThreshold(data, param = c("TTP_t0","MTT_t0"), threshold = 1:10,
                    hemisphere = "left", rm.CSF = TRUE,
                    GRalgo = TRUE, seed = c("MASK_T2_FLAIR_t2","MASK_DWI_t0"), W = W)

multiplot(selectCoords(MRIaggr.Pat1_red), main = "TTP_t0_thCC  - slice",
          hypoCC_Th1_10[,1], legend = FALSE, palette = rainbow(10), breaks = (0:10) - 10^{-10})

		}
}
\references{
Chantal Revol and Michel Jourlin. \emph{A new minimum varance region growing algorithm for image segmentation}. Pattern Recognition Letters, 18(3):249-258,1997.
}
