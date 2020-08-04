\docType{class}
\name{Carto3D-class}
\title{class "Carto3D"}
\alias{Carto3D-class}
\alias{Carto3D}
\description{
  Patient-specific storage of univariate contrast data.
}
\arguments{
  \item{identifier}{the patient identifier. \emph{character}.}  
  \item{parameter}{the name of the contrast parameter. \emph{character}.}   
  \item{contrast}{the contrast and the spatial coordinates of each voxel. \emph{data.frame} with four columns named \code{"i"} \code{"j"} \code{"k"} and the name of the contrast parameter.}  
  \item{fieldDim}{the dimension of the lattice containing the observations expressed in number of voxels. \emph{data.frame}.}  
  \item{voxelDim}{the voxel size with its unit of measure. A four columns \emph{data.frame} with names \code{"i"}, \code{"j"}, \code{"k"} and \code{"unit"}.}  
  \item{default_value}{the reference values of the contrast parameters (e.g. the background value). \emph{character}.}   
}

\section{S4 methods}{
Getters : 
\describe{
 \item{\code{\link{selectContrast}}}{return the contrast parameters.}
 \item{\code{\link{selectCoords}}}{return the coordinates.}
 \item{\code{\link{selectIdentifier}}}{return the identifier of the patient}
 \item{\code{\link{selectN}}}{return the number of voxels}
 \item{\code{\link{selectParameter}}}{return the names of the parameters stored in the object}
 \item{\code{\link{selectFieldDim}}}{return the dimension of the lattice containing the observations}
 \item{\code{\link{selectVoxelDim}}}{return the dimension of a voxel}
}

Displayers : 
\describe{
\item{\code{multiplot}}{display the contrast data.}
}

Constructors : 
\describe{
\item{\code{constMRIaggr}}{aggregate several Carto3D objects into a single objet MRIaggr.}
}

Initializers : 
\describe{
\item{\code{initNum}}{check the validity of the num argument.}
}

}
\seealso{
\code{\linkS4class{MRIaggr}} class that aggregates several \code{carto3D} objects.
\code{\link{Carto3D2MRIaggr}} to convert \code{carto3D} objects into \code{\linkS4class{MRIaggr}} objects.
}
\keyword{classes}
