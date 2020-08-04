\name{calcBlockW_cpp}
\title{Find disjoint spatial blocks of sites}
\alias{calcBlockW_cpp}
\description{
C++ function called by \code{\link{calcBlockW}} to identify independant spatial blocks of sites according to a neighbourhood matrix.
}
\usage{
calcBlockW_cpp(W_i, W_p, site_order, dist_center, dist_max, verbose)
}
\arguments{
  \item{W_i}{the 0-based row numbers for each non-zero element in the sparse neighbourhood  matrix. \emph{integer vector}.}  
  \item{W_p}{the pointers to the initial (zero-based) index of elements in the colum of the sparse neighbourhood matrix. \emph{integer vector}.}
  \item{site_order}{a specific order to go all over the sites. \emph{integer vector}.}  
  \item{dist_center}{the distance between each point and a reference point. \emph{numeric vector}.}  
  \item{dist_max}{the neighbourhood range. \emph{numeric vector}.}  
  \item{verbose}{Should the process be verbose over iterations ? \emph{logical}.}
}
\concept{cpp.}
\keyword{functions,internal}

