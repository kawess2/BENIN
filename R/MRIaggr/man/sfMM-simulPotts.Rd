\name{simulPotts}
\title{Draw a sample from a Potts model}
\alias{simulPotts}
\description{
Draw a sample from a Potts model. Call the \code{simulPotts_cpp} or the \code{simulPottsFast_cpp} C++ function.
}
\usage{
simulPotts(W, G, rho, initialization = NULL, site_order = NULL, iter_max = 2500,
         cv.criterion = 0, fast = FALSE, verbose = TRUE)
}
\arguments{
\item{W}{the neighbourhood matrix. \emph{dgCMatrix}. Should be normalized by row (i.e. \code{rowSums(W)=1}). REQUIRED.}
\item{G}{the number of groups. \emph{postive integer}. REQUIRED.}
\item{rho}{the value of the regularization parameter (also called potts parameter or inverse temperature). REQUIRED.}
\item{initialization}{the probability membership of each voxel to each group used for initialisation. \emph{matrix}.}
\item{site_order}{a specific order to go all over the sites. Can be \code{NULL}, \code{TRUE} or an \emph{integer vector}.}  
\item{iter_max}{the number of iterations. \emph{numeric}.}
\item{cv.criterion}{maximum difference in group probability membership for convergence ? \emph{double} between 0 and 1.}
\item{fast}{should \code{simulPottsFast_cpp} be used ? \emph{logical}. Faster but diseable tracing and convergence checking.}
\item{verbose}{should the simulation be be traced over iterations ? \emph{logical}.}
}
\details{
If no specific order is set, the order is sampled in a uniformed law which increase the computation time. 
If \code{site_order} is set to true, independant spatial block are identifyied using the \code{calcBlockW} function.
Each independant spatial block is then iteratively considered. 
}
\value{
A \emph{numeric vector} containing the regional potential.
}
\examples{
# spatial field
\dontshow{
n <- 10
iter_max <- 100
}
\dontrun{
n <- 30
iter_max <- 500
}
G <- 3
coords <- data.frame(which(matrix(0, nrow = n * G, ncol = n * G) == 0, arr.ind = TRUE), 1)
optionsMRIaggr(quantiles.legend = FALSE, axes = FALSE, num.main = FALSE)

# neighbourhood matrix
resW <- calcW(coords, range = sqrt(2), row.norm = TRUE, calcBlockW = TRUE)
W <- resW$W

# with no initial sample
res1 <- simulPotts(W, G = 3, rho = 3, iter_max = iter_max)
multiplot(coords,
          apply(res1$simulation, 1, which.max),
          breaks=seq(0.5, 3.5, 1))

res2 <- simulPotts(W, G = 4, rho = 3, iter_max = iter_max)
multiplot(coords,
          apply(res2$simulation, 1, which.max),
          breaks = seq(0.5, 4.5, 1))

res3 <- simulPotts(W, G = 4, rho = 6, iter_max = iter_max)
multiplot(coords,
          apply(res3$simulation, 1, which.max),
          breaks = seq(0.5, 4.5, 1))

# with specific initialisation
res3.bis <- simulPotts(W, rho = 6, initialization = res3$simulation, iter_max = iter_max / 2)
multiplot(coords,
          apply(res3.bis$simulation, 1, which.max),
          breaks=seq(0.5, 4.5, 1))

res3.ter <- simulPotts(W, rho = 6, initialization = res3$simulation, iter_max = iter_max)
multiplot(coords,
          apply(res3.ter$simulation, 1, which.max),
          breaks=seq(0.5, 4.5, 1))
		  
#### defining site order save time
site_order <- unlist(resW$blocks$ls_groups) - 1

system.time(
  res <- simulPotts(W, iter_max = 100, G = 3, rho = 6, site_order = NULL, 
                    fast = TRUE, verbose = FALSE)
)

system.time(
  res <- simulPotts(W, iter_max = 100, G = 3, rho = 6, site_order = site_order, 
                    fast = TRUE, verbose = FALSE)
)

system.time(
  res <- simulPotts(W, iter_max = 100, G = 3, rho = 6, site_order = NULL, 
                    fast = FALSE, verbose = FALSE)
)

system.time(
  res <- simulPotts(W, iter_max = 100, G = 3, rho = 6, site_order = site_order, 
                    fast = FALSE, verbose = FALSE)
)


}
