\name{constReduceMRIaggr}
\title{Reduce a MRIaggr}
\alias{constReduceMRIaggr}
\alias{constReduceMRIaggr,MRIaggr-method}
\description{
Construct a \code{\linkS4class{MRIaggr}} object restricted to a subset of observations.
}
\usage{
\S4method{constReduceMRIaggr}{MRIaggr}(object, mask, numeric2logical = FALSE, keep.index = TRUE)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{mask}{the binary contrast parameter or a vector indicating the observations to be kept. \emph{character} or \emph{logical vector} with length equal to the number of observations in \code{object}. REQUIRED.}
  \item{numeric2logical}{should \code{mask} be converted to logical ? \emph{logical}.}
  \item{keep.index}{should the previous \code{index} parameter be saved in the \code{ls_descStats} slot ? \emph{logical}.}
}
\value{
  a \code{\linkS4class{MRIaggr}} object.
}
\seealso{
\code{\link{calcBrainMask}} to compute an indicator of the brain observations.
}
\examples{
## load NIFTI files and convert them to MRIaggr
path.Pat1 <- system.file(file.path("nifti"), package = "MRIaggr")
ls.array <- list(readMRI(file.path(path.Pat1, "T2_GRE_t0"), format = "nifti"))
MRIaggr.Pat1 <- constMRIaggr(ls.array, identifier = "Pat1", param = "T2_GRE_t0")

## create the cerebral mask
\dontshow{
res <- calcBrainMask(MRIaggr.Pat1, param = "T2_GRE_t0", type = "threshold",
                     th.select_optima = 2, update.object = TRUE, overwrite = TRUE)
}
\dontrun{
res <- calcBrainMask(MRIaggr.Pat1, param = "T2_GRE_t0", type = "kmeans",
                     kmeans.n_groups = 2:4,
                     update.object = TRUE, overwrite = TRUE)

res <- calcSmoothMask(MRIaggr.Pat1, update.object = TRUE, overwrite = TRUE)
}
## display
multiplot(MRIaggr.Pat1, param = "mask", legend = FALSE)					  

## construct the reduced object
MRIaggr.Pat1_red <- constReduceMRIaggr(MRIaggr.Pat1, mask = "mask")

## display
par(mfrow = c(2,4), mar = rep(1.75, 4), mgp = c(2,0.75,0))
multiplot(MRIaggr.Pat1, param = "T2_GRE_t0",
             window = NULL,breaks = seq(0, 300, 1), legend = FALSE)
multiplot(MRIaggr.Pat1_red, param = "T2_GRE_t0",
             window = NULL, breaks = seq(0, 300, 1), legend = FALSE)
}
\concept{const.}
\keyword{methods}

