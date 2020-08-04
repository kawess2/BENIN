\name{calcPotts_cpp}
\title{Iterated conditional means for spatial regularization}
\alias{calcPotts_cpp}
\description{
C++ function that performs spatial regularization of probabilistic group membership using Iterated Conditional Means. 
For internal use.
}
\usage{
calcPotts_cpp(W_SR, W_LR, sample, rho, coords, site_order, iter_max, cv_criterion, 
         test_regional, distance_ref, threshold, neutre, nbGroup_min, multiV,
	     last_vs_others, prior, type_reg, verbose)
}
\arguments{
\item{W_SR}{the local neighborhood matrix. \emph{dgCMatrix}. Should be normalized by row (i.e. \code{rowSums(W_SR)=1}).}
\item{W_LR}{the regional neighborhood matrix. \emph{dgCMatrix}. Should contain the distances between the observations (\code{0} indicating infinite distance).}
\item{sample}{the initial group probability membership. \emph{numeric vector}.}
\item{rho}{the value of the spatial regularisation parameters. \emph{numeric vector}.}
\item{coords}{the voxel coordinates. \emph{matrix}.}
\item{site_order}{a specific order to go all over the sites. \emph{integer vector}.}  
\item{iter_max}{the maximum number of iterations. \emph{integer}.}
\item{cv_criterion}{the convergence criterion. \emph{numeric}.}

\item{test_regional}{should regional regularisation be used. \emph{logical}.}
\item{distance_ref}{the distance defining the several neighborhood orders relatively to \code{W_LR}. \emph{numeric vector}.}
\item{threshold}{the minimum value to consider non-negligible group membership. \emph{numeric}.}
\item{neutre}{the value to which the regional potential is set when the group is too small to compute it. \emph{numeric}.}
\item{nbGroup_min}{the minimum group size of the spatial groups required for performing regional regularization. \emph{integer}.}
\item{multiV}{should the regional potential range be specific to each spatial group ? \emph{logical}.}

\item{last_vs_others}{should the regional potential be computed for each class or only for the last versus the others ? \emph{logical}.}
\item{type_reg}{the type image to regularize : continuous or categorical ? \emph{logical}.}
\item{prior}{should the \code{sample} values be used as a prior ? \emph{logical}.}
\item{verbose}{should the running be be traced over iterations ? \emph{logical}. 1 to display each iteration and 2 to display convergence diagnostics.}
}
\seealso{
\code{\link{calcPotts}} which is the R interface to this C++ function.
}
\concept{cpp.}
\keyword{functions,internal}
              