\name{rhoMF}
\title{Estimation of the local and regional spatial correlation}
\alias{rhoMF}
\usage{
rhoMF(Y, W_SR, rho_max = 50, prior_prevalence = TRUE, 
    test.regional = FALSE, W_LR, distance.ref, coords, threshold = 0.1,
	nbGroup_min = 100, regionalGroups = "last_vs_others", multiV = TRUE)
}
\description{
Estimation the spatial regularization parameters on external data using mean field approximation.
}
\arguments{
  \item{Y}{a \emph{matrix} containing the observations (by rows) for the various groups (by columns). REQUIRED.}  
  \item{W_SR}{the local neighbourhood matrix. \code{dgCMatrix}. Should be normalized by row (i.e. \code{rowSums(Wweight_SR)=1}). REQUIRED.}
  \item{rho_max}{Maximum possible rho value (\emph{numeric}), minimum is 0.}
  \item{prior_prevalence}{should a prior on class prevalence be including when estimating the regularisation parameters ? \emph{logical}.}
  
  \item{test.regional}{Should regional regularization be considered. \emph{logical}.}
  \item{W_LR}{the regional neighbourhood matrix. \emph{dgCMatrix}. Should be contains the distances between the observations (0 indicating infinite distance).}
  \item{distance.ref}{the interval of distance defining the several neighbourhood orders in \code{W_LR}. \emph{numeric vector}.}
  \item{threshold}{the minimum value of the posterior probability for group G for being considered as lesioned when forming the spatial groups. \emph{double}.} 
  \item{nbGroup_min}{the minimum group size of the spatial groups required for computing the regional potential. \emph{integer}.}
  \item{coords}{coordinates of each site. \emph{matrix}.}
  \item{regionalGroups}{how should the regional potential be computed : last group versus the others (\code{"last_vs_others"}) or for each group (\code{"each"}).}
  \item{multiV}{should the regional neighbourhood range be computed for each spatial group ? \emph{logical}.}
}
\value{
  A \emph{numericVector} containing the estimated regularisation parameter(s).
}
\seealso{
\code{\link{calcW}} to compute the neighbourhood matrix, \cr
\code{\link{simulPotts}} to draw simulations from a Potts model. \cr
\code{\link{rhoLvfree}} to estimate the regularization parameters using mean field approximation.
\code{\link{calcPottsParameter}} general interface for estimating the regularization parameters.
}
\examples{
# spatial field
\dontrun{
n <- 50
}
\dontshow{
n <- 10
}
G <- 3
coords <- which(matrix(0, nrow = n * G, ncol = n * G) == 0,arr.ind = TRUE)

# neighbourhood matrix
W_SR <- calcW(as.data.frame(coords), range = sqrt(2), row.norm = TRUE)$W
W_LR <- calcW(as.data.frame(coords), range = 10, row.norm = FALSE)$W

# initialisation
set.seed(10)
sample <- simulPotts(W_SR, G = 3, rho = 3.5, iter_max = 500, 
                     site_order = TRUE)$simulation

multiplot(as.data.frame(coords), sample,palette = "rgb")

# estimation
rho <- rhoMF(Y=sample, W_SR = W_SR)
rho

# the regional potential is computed for each group
rho <- rhoMF(Y = sample, W_SR = W_SR,
             test.regional = TRUE, W_LR = W_LR, distance.ref = seq(1, 10, 0.5),
			 coords = coords, regionalGroups = "each")
rho

# the regional potential is computed for the last group vs. the others
rho <- rhoMF(Y = sample, W_SR = W_SR,
             test.regional = TRUE, W_LR = W_LR, distance.ref = seq(1, 10, 0.5),
			 coords = coords, regionalGroups = "last_vs_others")
rho
}



                        
                        
                      