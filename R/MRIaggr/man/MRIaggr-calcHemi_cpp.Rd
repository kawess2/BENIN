\name{calcHemi_cpp}
\title{Mid-saggital plan search}
\alias{calcHemi_cpp}
\description{
C++ function called by \code{\link{calcHemisphere}} to compute the objective function that evaluates the consistency of the mid-sagittal plan. For internal use.
}
\usage{
calcHemi_cpp(coordsI, coordsJ, ls_indexK, n_num, value, n,
         i_pos, j_pos, angle_pos,
         penaltyNA, sd_data, p, symetrie) 
}
\arguments{
  \item{coordsI}{the coordinates of each voxel along the first coordinate. \emph{integer vector}.}
  \item{coordsJ}{the coordinates of each voxel along the second coordinate. \emph{integer vector}.}
  \item{ls_indexK}{the index of the voxels on each slice. \emph{list}.}
  \item{n_num}{the number of slices. \emph{integer}.}
  \item{value}{the contrast parameter value for each voxel. \emph{numeric vector}.}
  \item{n}{the number of voxels. \emph{integer}.}
  \item{i_pos}{the first coordinate of the center of the mid-sagittal plane. \emph{numeric}.}
  \item{j_pos}{the second coordinate of the center of the mid-sagittal plane. \emph{numeric}.}
  \item{angle_pos}{the angle between the axis of the second coordinate and the mid-sagittal plane. \emph{numeric}.}
  \item{penaltyNA}{a penalty factor for the percent of voxel without controlateral correspondant. \emph{numeric}.}
  \item{sd_data}{the standard deviation of the parameter. \emph{numeric}.}
  \item{p}{the penalization factor. \emph{positive numeric}.}
  \item{symetrie}{the type of objective function. \code{TRUE} correspond to \code{"symmetry"} and \code{FALSE} to \code{"asymmetry"}.}
}
\concept{cpp.}
\keyword{functions,internal}
