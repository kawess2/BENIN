\name{simulPottsFast_cpp}
\title{Potts model simulation}
\alias{simulPottsFast_cpp}
\description{
Simulation of isotropic Potts model with local interaction by Gibbs sampling. Slightly faster compared to \code{\link{simulPotts_cpp}}.
}
\usage{
simulPottsFast_cpp(W_i, W_p, W_x, site_order, sample, rho, n, p, iter_nb)
}
\arguments{
  \item{W_i}{the 0-based row numbers for each non-zero element in the sparse neighbourhood  matrix. \emph{integer vector}.}  
  \item{W_p}{the pointers to the initial (zero-based) index of elements in the colum of the sparse neighbourhood matrix. \emph{integer vector}.}
  \item{W_x}{the non-zero elements of the matrix. \emph{numeric vector}.}  
  \item{site_order}{a specific order to go all over the sites. An \emph{integer vector} or a negative number to indicate no specific order.}
  \item{sample}{A \emph{matrix} containing the initial observations (by rows) for the various groups (by columns).}  
  \item{rho}{A \emph{numeric} corresponding to the Potts model parameter (i.e. intensity of the neighbourhood correlations).}
  \item{n}{The number of observations. \emph{integer}.}
  \item{p}{The number of groups or colors. \emph{integer}.}
  \item{iter_nb}{The maximum number of iteration of the Gibbs sampler. \emph{integer}.}
}
\concept{cpp.}
\keyword{functions,internal}


                        
                        
                      