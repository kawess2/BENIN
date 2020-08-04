\name{filtrage3Dmed_cpp}
\title{3D median filtering}
\alias{filtrage3Dmed_cpp}
\description{
C++ function call by \code{\link{calcFilter}} that performs three dimensional median filtering. For internal use.
}
\usage{
filtrage3Dmed_cpp (Vec_data, p_data, Vec_operateur, p_operateur, index_data, na_rm)
}
\arguments{
  \item{Vec_data}{vector of data to which the filter will be applied.}
  \item{p_data}{spatial dimensions of the data.}
  \item{Vec_operateur}{vector representing the filter to be applied.}
  \item{p_operateur}{spatial dimensions of the filter.}
  \item{index_data}{index of the non NA data.}
  \item{na_rm}{should the observations with missing values in their neighbourhood be removed ? Otherwise the ponderation is adjusted.}  
}
\concept{cpp.}
\keyword{functions,internal}