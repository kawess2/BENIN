\name{calcMultiPotential_cpp}
\title{Computation of the spatial potential}
\alias{calcMultiPotential_cpp}
\description{
C++ function that computes the regional potential of spatial groups.
For internal use.
}
\usage{
calcMultiPotential_cpp(W_SR, W_LR, sample, threshold, coords, 
         distance_ref, nbGroup_min, multiV, neutre)
}
\arguments{
\item{W_SR}{The local neighborhood matrix. \emph{dgCMatrix}. Should be normalized by row (i.e. \code{rowSums(W_SR)=1}).}
\item{W_LR}{The regional neighborhood matrix. \emph{dgCMatrix}. Should contain the distances between the observations (\code{0} indicating infinite distance).}
\item{sample}{The group probability membership. \emph{numeric vector}.}
\item{threshold}{The minimum value to consider non-negligible group membership. \emph{numeric}.}
\item{coords}{The voxel coordinates. \emph{matrix}.}
\item{distance_ref}{The intervals of distance defining the several neighborhood orders in \code{W_LR}. \emph{numeric vector}.}
\item{nbGroup_min}{The minimum group size of the spatial groups required for performing regional regularization. \emph{integer}.}
\item{multiV}{Should the regional potential range be specific to each spatial group ? \emph{logical}.}
\item{neutre}{the value to which the regional potential is set when the group is too small to compute it. \emph{numeric}.}
}
\seealso{
\code{\link{calcMultiPotential}} which is the R interface to this C++ function.
}
\concept{cpp.}
\keyword{functions,internal}
		