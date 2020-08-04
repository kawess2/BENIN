\name{calcContro_cpp}
\title{Compute contralateral normalization}
\alias{calcContro_cpp}
\description{
C++ function called by \code{\link{calcContralateral}} to seek the contralateral voxel and the contralateral values. For internal use.
}
\usage{
calcContro_cpp(contrast, coords_px, index_k, index_k_contro, d_lim, 
         lambda, param_ref, var_ref, type_moy, type_med, type_NN, verbose)
}
\arguments{
  \item{contrast}{the contrast value of each voxel on a given slice. \emph{matrix}.}
  \item{coords_px}{the coordinates of the observations. \emph{matrix}.}
  \item{index_k}{the index of the observations on the hemisphere of interest. \emph{integer vector}.}
  \item{index_k_contro}{the index of the observations on the contralateral hemisphere. \emph{integer vector}.}
  \item{d_lim}{the distance within which the contralateral values are considered. \emph{numeric}.}
  \item{lambda}{the importance of the penalization by the distance. \emph{numeric}.}
  \item{param_ref}{the parameter to be used as a reference to identify the contralateral voxel. \emph{character}.}
  \item{var_ref}{the variance of the reference parameter. \emph{numeric}.}
  \item{type_moy}{should the mean contralateral value be used ? \emph{logical}.}
  \item{type_med}{should the median contralateral value be used ? \emph{logical}.}
  \item{type_NN}{should the closest contralateral voxel according to the reference parameter be used ? \emph{logical}.}
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
}
\concept{cpp.}
\keyword{functions,internal}
