\name{filtrage2Dmed_cpp}
\title{2D median filtering}
\alias{filtrage2Dmed_cpp}
\description{
C++ function call by \code{\link{calcFilter}} that performs two dimensional median filtering. For internal use.
}
\usage{
filtrage2Dmed_cpp(M_data, M_operateur, index_data, na_rm)
}
\arguments{
  \item{M_data}{matrix to which the filter will be applied.}
  \item{M_operateur}{the filter to be applied.}
  \item{index_data}{index of the non NA data.}
  \item{na_rm}{should the observations with missing values in their neighbourhood be removed ? Otherwise the ponderation is adjusted.}  
}
\concept{cpp.}
\keyword{functions,internal}