\name{filtrage3D_cpp}
\title{3D filtering}
\alias{filtrage3D_cpp}
\description{
C++ function call by \code{\link{calcFilter}} that performs three dimensional filtering. For internal use.
}
\usage{
filtrage3D_cpp(Vec_data, p_data, Vec_operateur, p_operateur, index_data, 
         bilateral, na_rm)
}
\arguments{
  \item{Vec_data}{vector of data to which the filter will be applied.}
  \item{p_data}{spatial dimensions of the data.}
  \item{Vec_operateur}{vector representing the filter to be applied.}
  \item{p_operateur}{spatial dimensions of the filter.}
  \item{index_data}{index of the non NA data.}
  \item{bilateral}{should the influence of each neighbor be ponderated by the difference in signal with the considered observation ?}
  \item{na_rm}{should the observations with missing values in their neighbourhood be removed ? Otherwise the ponderation is adjusted.}  
}
\concept{cpp.}
\keyword{functions,internal}