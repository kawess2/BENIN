\name{calcGR}
\title{Interface to the Growing Region algorithm}
\alias{calcGR}
\description{
Call the \code{\link{GRalgo}} function to perform the Growing Region algorithm.
}
\usage{
calcGR(contrast, W, seed, sigma_max, range = c(-Inf, +Inf), range.seed = c(-Inf, +Inf),
         breaks = 100, rescale = FALSE, iter_max = 100, sd.robust = FALSE,
		 keep.lower = FALSE, keep.upper = FALSE, verbose = TRUE,
		 history.sigma = FALSE, history.step = FALSE, history.front = FALSE)
}
\arguments{
  \item{contrast}{the contrast value of each observation. \emph{numeric vector}. REQUIRED.}
  \item{W}{the neighbourhood matrix. \emph{dgCMatrix}. REQUIRED.}
  \item{seed}{the index of the initial seeds or a binary indicator of the initial seeds. \emph{positive integer vector} or \emph{logical vector}. REQUIRED.}
  \item{sigma_max}{the maximum admissible value for the variability of the group contrast. \emph{positive numeric}. REQUIRED.}
  \item{range}{the range of acceptable contrast values for the growing region group. \emph{numeric vector of size 2}.}
  \item{range.seed}{the range of acceptable contrast values for the seeds. \emph{numeric vector of size 2}.}
 
  \item{breaks}{the break points or the number of break points to use to categorize the contrast distribution. \emph{numeric vector} or \code{postive integer}.}
  \item{rescale}{should the contrast be scaled ? \emph{logical}.}  
  \item{iter_max}{the maximum number of iterations for the expansion of the growing region. \code{postive integer}.}  
  \item{sd.robust}{should the median absolute deviation be used to estimte the variability of the group contrast, or the standard deviation ? \emph{logical}.}
  \item{keep.lower}{should removing observations with high intensity of the region be forbidden ? \emph{logical}.}  
  \item{keep.upper}{should removing observations with low intensity of the region be forbidden ? \emph{logical}.}  
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}  
  
  \item{history.sigma}{should the values of sigma be recorded ? \emph{logical}.}
  \item{history.step}{should the number of observations included in the growing region set be recorded ? \emph{logical}.}
  \item{history.front}{should the propagation front of the GR set be recorded ? \emph{logical}.}
}
\details{
FUNCTION: \cr
This implementation of the Growing Region algorithm was been proposed by (Revol et al. 1997).
}
\value{
  An \emph{list} containing :
  \itemize{
  \item\code{[[GR]]} : the index of the observations in the growing region. \emph{integer vector}.
  \item\code{[[test.break]]} : whether the GR algorithm was interrupted an during execution. \emph{logical}.
  \item\code{[[iter]]} : the number of the last iteration of the algorithm. \emph{integer}.
  \item\code{[[test.id]]} : whether the GR set has stabilised during the last iteration. \emph{logical}.
  \item\code{[[sigma]]} : if \code{history.sigma} was set to \code{TRUE}, the value of the homogeneity criterion at the begining and the end of each step (in columns) for all steps (in row). \emph{numeric matrix}.
  \item\code{[[history_GR]]} : if \code{history.step} was set to \code{TRUE}, the step when each GR observation was included in the GR set. \emph{integer vector}.
  \item\code{[[breaks]]} : if \code{history.front} was set to \code{TRUE}, the values used to categorize the contrast. \emph{numeric vector}.
 }
}
\seealso{ 
 \code{\link{calcCriteriaGR}} for an automatic estimate of the sigma value.
}
\references{
Chantal Revol and Michel Jourlin. \emph{A new minimum varance region growing algorithm for image segmentation}. Pattern Recognition Letters, 18(3):249-258,1997.
}
\examples{
## load a MRIaggr object
data(MRIaggr.Pat1_red, package = "MRIaggr")

calcThresholdMRIaggr(MRIaggr.Pat1_red, param = c("TTP_t0","MTT_t0"), threshold = 1:10,
                     name_newparam = c("TTP.th_t0","MTT.th_t0"),
                     update.object = TRUE, overwrite = TRUE)
					 
## display raw parameter
multiplot(MRIaggr.Pat1_red, param = "TTP.th_t0", num = 3, numeric2logical = TRUE,
          index1 = list(coords = "MASK_DWI_t0", outline = TRUE))

## extract raw parameter, coordinates and compute the neighbourhood matrix
carto <- selectContrast(MRIaggr.Pat1_red, num = 3, hemisphere = "lesion",
                         param = c("TTP.th_t0","TTP_t0","MASK_DWI_t0"))
coords <- selectCoords(MRIaggr.Pat1_red, num = 3, hemisphere = "lesion")
W <- calcW(coords, range = sqrt(2))$W

## the seed is taken to be the point with the largest TTP in the lesion mask
indexN <- which(carto$MASK_DWI_t0 == 1)
seed <- indexN[which.max(carto[indexN,"TTP_t0"])] 

## Display step by step the GR algorithm with sigma = 1
for(iter in c(0,1,2,5,10)){
  resGR1 <- calcGR(contrast = carto$TTP.th_t0, W = W, 
                   seed = seed, sigma_max = 1, iter_max = iter, verbose = FALSE)
  
  multiplot(MRIaggr.Pat1_red, param = "TTP.th_t0", num = 3,hemisphere = "lesion", legend = FALSE,
         breaks = seq(0,10,0.1), numeric2logical = TRUE, cex = 2,
		 main = paste("iteration=",iter," - slice ", sep = ""),          
         index1 = list(coords = coords[resGR1$GR,], pch = 20, cex = 1),
         index2 = list(coords = coords[seed,], pch = 20, cex = 1)
  )
}

\dontrun{
## GR with sigma = 2.2
resGR2 <- calcGR(contrast = carto$TTP.th_t0, W = W, 
                 seed = seed, sigma_max = 2.2, iter_max = 50,
                 history.step = TRUE, history.front = TRUE)

## display 
# display the GR over the raw contrast
multiplot(MRIaggr.Pat1_red, param = "TTP.th_t0", num = 3, hemisphere = "lesion", legend = FALSE,
             breaks = seq(0,10,0.1), numeric2logical = TRUE, cex = 2,            
             index1 = list(coords = coords[resGR2$GR,], pch = 20, cex = 1)
)

# display the step of inclusion in GR group for each observation
multiplot(coords[resGR2$GR,],
             resGR2$history.step,breaks = 0:10,
             index1=list(coords = coords[seed,]),
             palette = rainbow(10)
)

# display the front propagation 
multiplot(coords[resGR2$GR,],
             resGR2$Mfront[,7],
             index1 = list(coords = coords[seed,])
)
}
}
\concept{calc.}
\keyword{functions}
