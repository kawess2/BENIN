\docType{class}
\name{MRIaggr-class}
\title{Class "MRIaggr"}
\alias{MRIaggr}
\alias{MRIaggr-class}
\description{
  Patient-specific storage of multivariate contrast data and clinical data.
}
\arguments{
  \item{identifier}{the patient identifier. \emph{character}.}  
  \item{contrast}{the value of the contrast parameters for each observation. \emph{data.frame}.}  
  \item{clinic}{the clinical data of the patient. \emph{data.frame}.}  
  \item{fieldDim}{the dimensions of the spatial field (expressed in number of voxels) containing the voxels. \emph{data.frame}.}  
  \item{voxelDim}{the voxel size dimensions its unit of measure. A four columns \emph{data.frame} with names \code{"i"}, \code{"j"}, \code{"k"} and \code{"unit"}.}  
  \item{default_value}{a reference value for each contrast parameter. \emph{data.frame}.}  

  \item{history}{a list of the \code{calc} or \code{const} methods that have been applied to the \code{MRIaggr} object . \emph{data.frame}.}  
  \item{normalization}{the normalization values of each contrast parameter. \emph{list}.}  
  \item{hemispheres}{the presence or absence of lesion(s) in each hemisphere. \emph{data.frame}.}  
  \item{midplane}{the position of the mid-sagittal plane. \emph{data.frame}.}  
  \item{W}{the neighbourhood matrix. \emph{list}.}  
 
  \item{table_lesion}{the vertical distribution of the lesion volumes. \emph{data.frame}.}  
  \item{table_reperfusion}{the volumic reperfusion data. \emph{data.frame}.}  
  \item{table_hypoperfusion}{the volumic hypoperfusion data. \emph{data.frame}.}  

  \item{ls_statDesc}{a slot to store additional data. \emph{list}.}  
}
\section{S4 methods}{
Getters : 
\describe{
 \item{\code{\link{selectContrast}}}{return the contrast parameters.}
 \item{\code{\link{selectCoords}}}{return the coordinates.}
 \item{\code{\link{selectClinic}}}{return the clinical data.}
 \item{\code{\link{selectDefault_value}}}{return the default values of the contrast parameters.}
 \item{\code{\link{selectDescStats}}}{return extra information.}
 \item{\code{\link{selectHemispheres}}}{return the position of the lesion in each hemisphere.}
 \item{\code{\link{selectHistory}}}{return the call of the methods applied on the object.}
 \item{\code{\link{selectIdentifier}}}{return the identifier of the patient}
 \item{\code{\link{selectMidplane}}}{return the position of the mid-sagittal plan}
 \item{\code{\link{selectN}}}{return the number of voxels}
 \item{\code{\link{selectNormalization}}}{return the normalization values of the contrast parameters}
 \item{\code{\link{selectParameter}}}{return the names of the parameters stored in the object}
 \item{\code{\link{selectTable}}}{return volumic information about the lesion or about the perfusion}
 \item{\code{\link{selectFieldDim}}}{return the dimension of the lattice containing the observations}
 \item{\code{\link{selectVoxelDim}}}{return the dimension of a voxel}
}

Setters : 
\describe{
\item{\code{\link{allocContrast<-}}}{allocate a contrast parameter}
\item{\code{\link{allocClinic<-}}}{allocate clinical data}
\item{\code{\link{allocDescStats<-}}}{allocate extra information}
\item{\code{\link{allocHemisphere<-}}}{allocate the position of the lesion in the hemisphere or the position of the mid-sagittal brain}
\item{\code{\link{allocNormalization<-}}}{allocate normalization values}
\item{\code{\link{allocTable<-}}}{allocate volumic information}
\item{\code{\link{supprContrast<-}}}{remove a contrast parameter}
\item{\code{\link{supprDescStats<-}}}{remove extra information}
}

Calculators : 
\describe{
\item{\code{\link{calcBrainMask}}}{}
\item{\code{\link{calcContralateral}}}{}
\item{\code{\link{calcDistMask}}}{}
\item{\code{\link{calcDistTissues}}}{}
\item{\code{\link{calcFilter}}}{}
\item{\code{\link{calcGroupsMask}}}{}
\item{\code{\link{calcHemisphere}}}{}
\item{\code{\link{calcROCthreshold}}}{}
\item{\code{\link{calcNormalization}}}{}
\item{\code{\link{calcRegionalContrast}}}{}
\item{\code{\link{calcSmoothMask}}}{}
\item{\code{\link{calcTableHypoReperf}}}{}
\item{\code{\link{calcTableLesion}}}{}
\item{\code{\link{calcThresholdMRIaggr}}}{}
\item{\code{\link{calcTissueType}}}{}
\item{\code{\link{calcW}}}{}
}

Displayers : 
\describe{
\item{\code{\link{boxplotMask}}}{}
\item{\code{\link{heatmapMRIaggr}}}{}
\item{\code{\link{multiplot}}}{}
\item{\code{\link{plotDistClass}}}{}
\item{\code{\link{pointsHemisphere}}}{}
\item{\code{\link{plotLesion3D}}}{}
\item{\code{\link{plotTableLesion}}}{}
\item{\code{\link{outlineMRIaggr}}}{}
\item{\code{\link{summary,MRIaggr-method}}}{}
}

Constructors : 
\describe{
\item{\code{\link{constCompressMRIaggr}}}{generate from an existing MRIaggr object a new one with lower spatial resolution}
\item{\code{\link{constReduceMRIaggr}}}{copy an existing MRIaggr object restricted to a subset of voxels}
\item{\code{\link{writeMRIaggr}}}{create an image file containing the values of a contrast parameter stored in the MRIaggr object.}
}

Initializers : 
\describe{
\item{\code{\link{initNum}}}{initialization of the num argument} 
\item{\code{\link{initParameter}}}{intialization of the param argument}  
\item{\code{\link{initMask}}}{initialization of the mask argument}  
}
}
\seealso{
\code{\link{constMRIaggr}} to build a \code{\linkS4class{MRIaggr}} object from a list of array.  \cr
\code{\link{readMRI}} to read image files. 
}
\keyword{classes}
