\name{calcSigmaGR}
\title{Automatic Growing Region algorithm}
\alias{calcSigmaGR}
\description{
Evaluate the quality of the Growing Region partition regarding several homogeneity parameters.
}
\usage{
calcSigmaGR(contrast, W, seed, sigma, criterion.transition = FALSE,
         criterion.sdfront = FALSE, criterion.entropy = TRUE, criterion.Kalinsky = TRUE,
         criterion.Laboure = TRUE, verbose = TRUE, ...)
}
\arguments{
  \item{contrast}{the contrast value of each observation. \emph{numeric vector}. REQUIRED.}
  \item{W}{the neighbourhood matrix. \emph{dgCMatrix}. REQUIRED.}
  \item{seed}{the index of the initial seeds or a binary indicator of the initial seeds. \emph{positive integer vector} or \emph{logical vector}. REQUIRED.}
  \item{sigma}{the sequence of maximum admissible values for the group variability \emph{positive numeric vector}. REQUIRED.}
  
  \item{criterion.transition}{should the boundary criterion based on the transition levels be computed ? \emph{logical}. }  
  \item{criterion.sdfront}{should the boundary criterion based on the standard deviation be computed ? \emph{logical}.}  
  \item{criterion.entropy}{should the region criterion based on the entropy be computed ? \emph{logical}.}  
  \item{criterion.Kalinsky}{should the region criterion based on the Kalinsky index be computed ? \emph{logical}.}  
  \item{criterion.Laboure}{should the region criterion based on the Laboure index  be computed ? \emph{logical}.}
  
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}  
  \item{...}{arguments to be passed to \code{\link{calcGR}} for specifying the settings of the growing region algorithm : \code{range}, \code{range.seed}, \code{breaks}, \code{scale}, \code{iter_max}, \code{sd.robust}, \code{keep.lower} and  \code{keep.upper}.}  

 }
\details{
ARGUMENTS: \cr
Information about the \code{window}, \code{filename}, \code{width}, \code{height}, \code{path}, \code{unit} and \code{res} arguments can be found in the details section of \code{\link{initWindow}}.

Information about the \code{mar} and \code{mgp} arguments can be found in \code{\link{par}}.

FUNCTION: \cr
This implementation of the automated Growing Region algorithm was proposed by (Revol et al. 2002) : \code{criterion.transition} corresponds to w2, \code{criterion.sdfront} corresponds to w3 where m=f(x),  \code{criterion.entropy} corresponds to S(sigma_{max}) and \code{criterion.Laboure} corresponds to InvDGL(sigma).
\code{criterion.Kalinsky} corresponds to the Kalinsky criterion which is the ratio of the variance between groups over the variance withing groups.
}
\value{
  An \emph{list} containing :
  \itemize{
  \item\code{[[df.criterion]]} : the value of the clustering criterion (in columns) for each sigma value (in rows). \emph{numeric matrix}.
  \item\code{[[list.GR]]} : the list of the optimal GR sets, one for each clustering criterion.
  \item\code{[[best]]} : the optimal value of each clustering criterion. \emph{data.frame}.
  \item\code{[[n.max]]} : the number of observations. \emph{integer}.
 }
}
\seealso{ 
 \code{\link{plotSigmaGR}} for a graphical display of the quality criteria.
}
\references{
Chantal Revol-Muller, Francoise Peyrin, Yannick Carrillon and Christophe Odet. \emph{Automated 3D region growing algorithm based on an assessment function}. Pattern Recognition Letters, 23:137-150,2002.
}
\examples{
\dontrun{
## load an \code{MRIaggr} object
data(MRIaggr.Pat1_red, package = "MRIaggr")

calcThresholdMRIaggr(MRIaggr.Pat1_red, param = c("TTP_t0","MTT_t0"), threshold = 1:10,
                     name_newparam = c("TTP.th_t0","MTT.th_t0"),
                     update.object = TRUE, overwrite = TRUE)

## display raw parameter
multiplot(MRIaggr.Pat1_red, param="TTP.th_t0", num = 3, numeric2logical = TRUE,
          index1 = list(coords = "MASK_DWI_t0", outline = TRUE))

## extract raw parameter, coordinates and compute the neighbourhood matrix
carto <- selectContrast(MRIaggr.Pat1_red, num = 3, hemisphere = "lesion",
                        param = c("TTP.th_t0","TTP_t0","MASK_DWI_t0"))
coords <- selectCoords(MRIaggr.Pat1_red, num = 3, hemisphere = "lesion")
W <- calcW(coords, range = sqrt(2))$W

## the seed is taken to be the point with the largest TTP in the lesion mask
indexN <- which(carto$MASK_DWI_t0 == 1)
seed <- indexN[which.max(carto[indexN,"TTP_t0"])]

## find optimal sigma
resGR_sigma <- calcSigmaGR(contrast = carto$TTP.th_t0, W = W, seed = seed,
                           sigma = seq(1,4,0.1), iter_max = 50,
                           keep.upper = TRUE)

## display quality criteria according to sigma
plotSigmaGR(resGR_sigma)

## display retained region
multiplot(MRIaggr.Pat1_red, param = "TTP.th_t0", num = 3, numeric2logical = TRUE,
          index1 = list(coords = coords[resGR_sigma$list.GR$entropy,], outline = TRUE))

multiplot(MRIaggr.Pat1_red, param = "TTP.th_t0", num = 3, numeric2logical = TRUE,
          index1 = list(coords = coords[resGR_sigma$list.GR$Kalinsky,], outline = TRUE))

## find optimal sigma
resGR_sigma <- calcSigmaGR(contrast = carto$TTP.th_t0, W = W, seed = seed,
                           sigma = seq(1,4,0.1), iter_max = 50,
                           keep.upper = TRUE, keep.lower = TRUE)

## display quality criteria according to sigma
plotSigmaGR(resGR_sigma)

## display retained region
multiplot(MRIaggr.Pat1_red, param = "TTP.th_t0", num = 3, numeric2logical = TRUE,
          index1 = list(coords = coords[resGR_sigma$list.GR$entropy,], outline = TRUE))

multiplot(MRIaggr.Pat1_red, param = "TTP.th_t0", num = 3, numeric2logical = TRUE,
          index1 = list(coords = coords[resGR_sigma$list.GR$Kalinsky,], outline = TRUE))
}
}
\concept{calc.}
\keyword{functions}
