\name{calcMultiPotential}
\title{Computation of the spatial potential}
\alias{calcMultiPotential}
\description{
Compute the regional potential of a spatial field. Call the \code{calcMultiPotential_cpp} C++ function.
}
\usage{
calcMultiPotential(W_SR, W_LR, distance.ref, sample, coords,
         threshold = 0.01, nbGroup_min = 100, check.args = TRUE, verbose = TRUE)
}
\arguments{
\item{W_SR}{The local neighborhood matrix. \emph{dgCMatrix}. Should be normalized by row (i.e. \code{rowSums(Wweight_SR)=1}). REQUIRED.}
\item{W_LR}{The regional neighborhood matrix. \emph{dgCMatrix}. Should contain the distances between the observations (\code{0} indicating infinite distance). REQUIRED.}
\item{distance.ref}{The intervals of distance defining the several neighborhood orders in \code{W_LR}. \emph{numeric vector}. REQUIRED.}
\item{sample}{The probability membership to the group. \emph{numeric vector} or \emph{character vector}. REQUIRED.}
\item{coords}{The voxel coordinates. \emph{matrix}. REQUIRED.}
\item{threshold}{The minimum value to consider non-negligible group membership. \emph{numeric}. Default is  \code{0.01}.}
\item{nbGroup_min}{The minimum group size of the spatial groups required for computing the potential. \emph{integer}. Default is \code{100}.}
\item{check.args}{Should the validity of the arguments be checked ? \emph{logical}.}
\item{verbose}{Should the radius of the spatial groups be displayed ? \emph{logical}.}
}
\details{
If \code{check.args} is set to \code{TRUE}, argument \code{ccords} must be a \emph{matrix} (and not a \emph{data.frame}) and the x slot of argument \code{W_LR} must contain categories of distances instead of distances.
Distance categories must begin at \code{0} and end at \code{length(distance.ref)-1}.
}
\value{
A \emph{numeric vector} containing the regional potential.
}
\examples{
# spatial field
\dontrun{
n <- 30
}
\dontshow{
n <- 10
}
G <- 3
coords <- data.frame(which(matrix(0, nrow = n * G, ncol = n * G) == 0, arr.ind = TRUE), 1)
optionsMRIaggr(quantiles.legend = FALSE, axes = FALSE, num.main = FALSE)

# neighborhood matrix
W_SR <- calcW(coords, range = sqrt(2), row.norm = TRUE)$W
W_LR <- calcW(coords, range = 10, row.norm = FALSE)$W
distance.ref <- seq(1, 10, 1)

# data
set.seed(10)
sample <- simulPotts(W = W_SR, G = G, rho = 3.5, iter_max = 500)

# 
multiplot(coords,sample$simulation[,2])

V1 <- calcMultiPotential(W_SR = W_SR, W_LR = W_LR, distance.ref = seq(1, 10, 1),
     sample = sample$simulation[,2], coords = coords, verbose = TRUE)

multiplot(coords, V1)		


#
sampleV <- (apply(sample$simulation, 1, which.max) - 1) / 2
multiplot(coords, sampleV)

system.time(
  V2 <- calcMultiPotential(W_SR = W_SR, W_LR = W_LR, distance.ref = seq(1, 10, 1),
     sample = sampleV, coords = coords, verbose = TRUE)
)

Wcat_LR <- W_LR
Wcat_LR@x <- findInterval(x = Wcat_LR@x, vec = distance.ref) - 1

system.time(
  V2 <- calcMultiPotential(W_SR = W_SR, W_LR = Wcat_LR, distance.ref = seq(1, 10, 1),
     sample = sampleV, coords = as.matrix(coords), verbose = TRUE, check.args = FALSE)
)
# quicker but arguments have to be correctly specified

multiplot(coords, V2)		
	
}
\keyword{functions}