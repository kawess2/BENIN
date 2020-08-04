\name{calcPotts}
\title{Spatial regularization using ICM}
\alias{calcPotts}
\description{
Interface to C++ functions that perform spatial regularization of probabilistic group membership using Iterated Conditional Means. 
}
\usage{
calcPotts(W_SR, sample, rho, prior = TRUE, site_order = NULL, W_LR = NULL, 
         nbGroup_min = 100, coords = NULL, distance.ref = NULL, threshold = 0.1,
		 multiV = TRUE, iter_max = 200, cv.criterion = 0.005, verbose = 2)
}
\arguments{
\item{W_SR}{The local neighborhood matrix. \emph{dgCMatrix}. Should be normalized by row (i.e. \code{rowSums(W_SR)=1}). REQUIRED.}
\item{sample}{The initial group probability membership. \emph{numeric vector}. REQUIRED.}
\item{rho}{Value of the spatial regularisation parameters. \emph{numeric vector}. REQUIRED.}
\item{prior}{Should the \code{sample} values be used as a prior ? \emph{logical}.}
\item{site_order}{a specific order to go all over the sites. \emph{integer vector}.}  

\item{W_LR}{The regional neighborhood matrix. \emph{dgCMatrix}. Should contain the distances between the observations (\code{0} indicating infinite distance).}
\item{nbGroup_min}{The minimum group size of the spatial groups required for performing regional regularization. \emph{integer}.}
\item{coords}{The voxel coordinates. \emph{matrix}.}
\item{distance.ref}{The intervals of distance defining the several neighborhood orders in \code{W_LR}. \emph{numeric vector}.}
\item{threshold}{The minimum value to consider non-negligible group membership. \emph{numeric}.}
\item{multiV}{Should the regional potential range be specific to each spatial group ? \emph{logical}.}

\item{iter_max}{Maximum number of ICM iterations. \emph{integer}.}
\item{cv.criterion}{Convergence criterion of the ICM . \emph{numeric}.}
\item{verbose}{should the ICM be be traced over iterations ? \emph{logical}. 1 to display each iteration and 2 to display convergence diagnostics.}
}
\details{
The convergence criterion of the ICM is computed as maximum absolute difference between the group membership probability between two consecutive iterations.
}
\examples{
optionsMRIaggr(legend=FALSE,axes=FALSE,num.main=FALSE,mar=c(0,2,2,0))

# spatial field
\dontrun{
n <- 30
}
\dontshow{
n <- 10
}
G <- 3
coords <- which(matrix(0, nrow = n * G, ncol = n * G) == 0,arr.ind = TRUE)

# neighborhood matrix
W_SR <- calcW(as.data.frame(coords), range = sqrt(2), row.norm = TRUE)$W
resW <- calcW(as.data.frame(coords), range = 10, row.norm = FALSE, calcBlockW = TRUE)
W_LR <- resW$W
site_order <- unlist(resW$blocks$ls_groups) - 1

# initialisation
set.seed(10)
system.time(
sample <- simulPotts(W_SR, G = 3, rho = 3.5, iter_max = 500, 
     site_order = site_order)$simulation
)

intensity <- rnorm((n * G)^2, mean = apply(sample, 1, which.max), sd = 0.5)
likelihood <- matrix(unlist(lapply(1:3, function(x){dnorm(intensity, mean = x, sd = 0.5)})),
     ncol = G, nrow = (n * G)^2, byrow = FALSE)
likelihood_sqrt <- sqrt(likelihood)

probability <- sweep(likelihood_sqrt, MARGIN = 1, FUN = "/", STATS = rowSums(likelihood_sqrt))
  
multiplot(as.data.frame(coords), probability, palette = "rgb",
     main = "original image")

#### local image restoration
LocalRestoration <- calcPotts(W_SR = W_SR, sample = probability, rho = 4,
     site_order = site_order)

multiplot(as.data.frame(coords), LocalRestoration$predicted, palette = "rgb",
     main = "local restoration of the image")

#### regional image restoration
distance.ref <- seq(1, 10, 1)

RegionalRestoration <- calcPotts(W_SR = W_SR, sample = probability, 
     rho = c(4,2), site_order = site_order,
     W_LR = W_LR, coords = coords, distance.ref = distance.ref)
  
# regional potential  
multiplot(as.data.frame(coords),
     matrix(unlist(RegionalRestoration$Vregional), ncol = 3, nrow = (n * G)^2, byrow = FALSE),
     palette = "rgb", main = "regional potentials")

# final image
multiplot(as.data.frame(coords),
     matrix(unlist(RegionalRestoration$predicted), ncol = 3, nrow = (n * G)^2, byrow = FALSE),
     palette = "rgb", main = "local and regional \n restoration of the image")
}

