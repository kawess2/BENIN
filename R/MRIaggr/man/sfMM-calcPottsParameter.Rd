\name{calcPottsParameter}
\title{Estimation of the local regularization parameters}
\alias{calcPottsParameter}
\usage{
calcPottsParameter(Y, W_SR, coords = NULL, range = NULL, method = "MF", verbose = 3, ...)
}
\description{
Estimation of the local regularization parameters using mean field approximation or a likelihood free method.
}
\arguments{
  \item{Y}{a \emph{matrix} containing the observations (by rows) for the various groups (by columns). REQUIRED.}  
  \item{W_SR}{the local neighbourhood matrix. \emph{dgCMatrix}. Should be normalized by row (i.e. \code{rowSums(W_SR)=1}). REQUIRED.}
  \item{coords}{if argument \code{W_SR} is not specified, the coordinates of the sites used to compute the local neighbourhood matrix. \emph{data.frame}.}
  \item{range}{if argument \code{W_SR} is not specified, the range of the neighbourhood used to compute the local neighbourhood matrix. positive \emph{double}.}
  \item{method}{the method used to estimate the regularization parameter. Can be either \code{"MF"} or \code{"Lvfree"}.}
  \item{verbose}{should the execution of the function be be traced ? Can be \code{0}, \code{1}, \code{2} or \code{3}.}  
  \item{...}{additional arguments to be passed to the function called by \code{calcPottsParameter} to perform the estimation. See the arguments of \code{link{rhoMF}} or \code{link{rhoLvfree}}.}
  }
\details{
FUNCTION: \cr
This function call either \code{rhoMF} if argument \code{method} is \code{"MF"} or \code{\link{rhoLvfree}} if argument \code{method} is \code{"Lvfree"}.
Estimation of the regional regularization parameter is only available with mean field approximation (\code{method="MF"}).
The likelihood free estimation (Pereyra et al, 2013) should give a more accurate estimation despite a higher computational cost.
}
\references{
M. Pereyra, N. Dobigeon, H. Batatia, and J.Y. Tourneret. \emph{Estimation the granularity coefficient of a Potts-Markov random field within an MCMC algorithm}. IEEE Trans. Image Porcessing, 22(6):2385-2397, 2013.
}
\seealso{
\code{\link{calcW}} to compute the neighbourhood matrix, \cr
\code{\link{simulPotts}} to simulate from a Potts model.
}
\examples{
# spatial field
\dontrun{
n <- 50
}
\dontshow{
n <- 5
}
G <- 3
coords <- as.data.frame(which(matrix(0, nrow = n * G, ncol = n * G) == 0, arr.ind = TRUE))

# neighbourhood matrix
resW <- calcW(as.data.frame(coords), range = sqrt(2), row.norm = TRUE, calcBlockW = TRUE)
W_SR <- resW$W
site_order <- unlist(resW$blocks$ls_groups) - 1

# initialisation
set.seed(10)
sample <- simulPotts(W_SR, G = 3, rho = 3.5, iter_max = 500, 
                     site_order = site_order)$simulation

#### estimation using the neighbourhood matrix
rho_MF <- calcPottsParameter(Y = sample, W_SR = W_SR,  method = "MF")

\dontrun{
rho_Lvfree <- calcPottsParameter(Y = sample, W_SR = W_SR, 
                                site_order = site_order, method = "Lvfree", verbose = 2)
}

#### estimation using the coordinates
rho_MF <- calcPottsParameter(Y = sample, coords = coords, range = sqrt(2), method = "MF")

\dontrun{
rho_Lvfree <- calcPottsParameter(Y = sample, coords = coords, range = sqrt(2),
                                site_order = site_order, method = "Lvfree", verbose = 2)
}
\dontshow{
rho_Lvfree <- calcPottsParameter(Y = sample, coords = coords, range = sqrt(2), n.sample = 500,
                                site_order = site_order, method = "Lvfree", verbose = 2)
}
}


                        
                        
                      
