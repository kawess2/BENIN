\name{calcGroupsW_cpp}
\title{Find spatial groups}
\alias{calcGroupsW_cpp}
\description{
C++ function called by \code{\link{calcGroupsW}} to compute the spatial groups. For internal use.
}
\usage{
calcGroupsW_cpp(W_i, W_p, subset, max_groups, verbose)
}
\arguments{
  \item{W_i}{the 0-based row numbers for each non-zero element in the sparse neighbourhood  matrix. \emph{integer vector}.}  
  \item{W_p}{the pointers to the initial (zero-based) index of elements in the colum of the sparse neighbourhood matrix. \emph{integer vector}.}
  \item{subset}{the subset of observations to use. \emph{positive integer vector}.}
  \item{max_groups}{the maximum number of groups. \emph{postive integer}.}
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
}
\concept{cpp.}
\keyword{functions,internal}


