\name{rhoLvfree}
\title{Estimation of the local regularization parameters}
\alias{rhoLvfree}
\usage{
rhoLvfree(Y, W_SR, rho_max = 10, rho_init = "PML", sd_rho_init = "PML", site_order = NULL,
     dprior = function(x){stats::dunif(x,0,rho_max)}, epsilon = 0.001, update_epsilon = 1,
	 n.sample = 2500, iter_max = 3, burn_in = round(0.25 * n.sample), thin = 1,
	 n.chain = 1, verbose = 3, cpus = 1, export.coda = TRUE)
}
\description{
Estimation of the local regularization parameters on external data using an efficient MCMC procedure.
}
\arguments{
  \item{Y}{a \emph{matrix} containing the observations (by rows) for the various groups (by columns). REQUIRED.}  
  \item{W_SR}{the local neighbourhood matrix. \emph{dgCMatrix}. Should be normalized by row (i.e. \code{rowSums(W_SR)=1}). REQUIRED.}
  \item{rho_max}{maximum possible rho value, minimum is 0. \emph{double}.}
  \item{rho_init}{the initial rho value. \emph{double}. If \code{NULL}, it is sample in a uniform law between 0 and \code{rho_max}}
  \item{sd_rho_init}{the standard deviation of the rho value for the proposal distribution. \emph{double}. It is updated at the end of the burn in.}
  \item{site_order}{a specific order to go all over the sites. Can be \code{NULL} or an \emph{integer vector}.}  
  
  \item{dprior}{the prior distribution. \emph{function}. Must be a density function between 0 and \code{rho_max}.}
  \item{epsilon}{the tolerance parameter regulating the acceptance of the parameter. \emph{double}.}
  \item{update_epsilon}{if over 1 the value by which the epsilon parameter is devided every 5\% of the total iterations, if under one the acceptance rate to reach. \emph{double}.}
  \item{n.sample}{the number of iterations of the Gibbs sampler. \emph{integer}.}
  \item{iter_max}{the number of Gibbs move for each simulation. \emph{integer}.}
  \item{burn_in}{the number of iteration of the burn in phase. \emph{integer}.}
  \item{thin}{the thinning interval between consecutive observations. \emph{integer}.}
  \item{n.chain}{the number of chain to use. \emph{integer}.}
  \item{verbose}{how should the execution of the function should be traced ? 1 traces display the initialization and the final result, \code{2} traces every 5\% of the total iterations and \code{3} also displays the update performed every 5\%.}
  \item{cpus}{the number of CPU to use. \emph{integer}.}
  \item{export.coda}{should the results be convert to the mcmc format of the coda package ? \emph{logical}.}
}
\details{
    FUNCTION: \cr
This function uses a likelihood-free Metropolis-Hastings method (proposed by Pereyra et al., 2013) to estimate the local regularization parameter.
It should give a less biased estimation of the parameter but require a higher computational cost.
}
\value{
   A good initialisation (\code{rho_init}) may be obtained with the \code{rhoMF} function. 
   It is automatically done if rho_init is set to \code{"PML"}. 
   
   Setting \code{sd_rho_init} to \code{"PML"} leads to adjustement of the initial variance according to the rho value and the sample size. 
}
\references{
M. Pereyra, N. Dobigeon, H. Batatia, and J.Y. Tourneret. \emph{Estimation the granularity coefficient of a Potts-Markov random field within an MCMC algorithm}. IEEE Trans. Image Porcessing, 22(6):2385-2397, 2013.
}
\seealso{
\code{\link{calcW}} to compute the neighbourhood matrix, \cr
\code{\link{simulPotts}} to simulate from a Potts model.
\code{\link{rhoMF}} to estimate the regularization parameters using mean field approximation.
\code{\link{calcPottsParameter}} general interface for estimating the regularization parameters.
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
coords <- which(matrix(0, nrow = n * G, ncol = n * G) == 0, arr.ind = TRUE)

# neighbourhood matrix
resW <- calcW(as.data.frame(coords), range = sqrt(2), row.norm = TRUE, calcBlockW = TRUE)
W_SR <- resW$W
site_order <- unlist(resW$blocks$ls_groups)-1

# initialisation
set.seed(10)
sample <- simulPotts(W_SR, G = 3, rho = 3.5, iter_max = 500, 
                     site_order = site_order)$simulation

multiplot(as.data.frame(coords), sample,palette = "rgb")

# estimation

\dontrun{
rho <- rhoLvfree(Y = sample, W_SR = W_SR, site_order = site_order)
}
\dontshow{
rho <- rhoLvfree(Y = sample, W_SR = W_SR, site_order = site_order, n.sample = 500)
}
}


                        
                        
                      
