\name{constMRIaggr}
\title{Array constructor for MRIaggr object}
\alias{constMRIaggr}
\description{
Construct a \code{\linkS4class{MRIaggr}} object from a list of array, each array corresponding to a different contrast parameters.
}
\usage{
constMRIaggr(ls.array, identifier, param, default_value = NULL,
     pos_default_value = c(1,1,1), tol = 10^{-10},
	 voxelDim = NULL, verbose = optionsMRIaggr("verbose"), rm.ls.array = FALSE)
}
\arguments{
  \item{ls.array}{the value of the contrast parameter(s) for each observation. \code{list of array}. REQUIRED.}
  \item{identifier}{the identifier of the patient to which belong the contrast parameters. \emph{character}. REQUIRED.}
  \item{param}{the contrast parameter(s). \emph{character vector} or \code{NULL}. REQUIRED.}
  \item{default_value}{the reference values of the contrast parameters (e.g. the background value). \emph{character} or \code{NULL} leading to search the reference value in \code{array[pos_default_value]}.}
  \item{pos_default_value}{the coordinates of the observations that contains the reference value. \emph{numeric vector}.}   
  \item{tol}{numeric precision for the consistency check. \emph{positive numeric}.}
  
  \item{voxelDim}{the voxel size with its unit of measure. A four columns \emph{data.frame} with names \code{"i"}, \code{"j"}, \code{"k"} and \code{"unit"}.}
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
  \item{rm.ls.array}{should the object on which the \code{ls.array} argument points be removed form the global environment ? \emph{logical}.}
}
\details{
ARGUMENTS: \cr
All the array in \code{ls.array} in must have the same dimensions.

Information about the \code{param} argument can be found in the details section of \code{\link{initParameter}}.

\code{pos_default_value} is active only if \code{default_value} is set to \code{NULL}.
}
\value{
  a \code{\linkS4class{MRIaggr}} object.
}
\examples{
#### 1- 1st method ####
## load NIFTI files
path.Pat1 <- system.file(file.path("nifti"), package = "MRIaggr")

nifti.Pat1_TTP_t0 <- readMRI(file.path(path.Pat1, "TTP_t0"), format = "nifti")
nifti.Pat1_DWI_t0 <- readMRI(file.path(path.Pat1, "DWI_t0"), format = "nifti")
nifti.Pat1_MASK_DWI_t0 <- readMRI(file.path(path.Pat1, "MASK_DWI_t0"), format = "nifti")
nifti.Pat1_MASK_T2_FLAIR_t2 <- readMRI(file.path(path.Pat1, "MASK_T2_FLAIR_t2"), 
                                       format = "nifti")

## convert them to MRIaggr
MRIaggr.Pat1 <- constMRIaggr(list(nifti.Pat1_TTP_t0, nifti.Pat1_DWI_t0,
                                  nifti.Pat1_MASK_DWI_t0, nifti.Pat1_MASK_T2_FLAIR_t2),
                identifier= "Pat1", param=c("TTP_t0","DWI_t0","MASK_DWI_t0","MASK_T2_FLAIR_t2"))


#### 2- 2nd method ####
## load nifti files
param <- c("DWI_t0.nii","MASK_DWI_t0.nii","MTT_t0.nii","TTP_t0.nii","T1_t0.nii","T2_GRE_t0.nii",
           "MTT_t1.nii","TTP_t1.nii","T2_FLAIR_t2.nii","MASK_T2_FLAIR_t2.nii")

ls.array <- list()
for(iter_param in 1:length(param)){
  ls.array[[iter_param]] <- readMRI(file.path(path.Pat1,param[iter_param]), format = "nifti")
}

## convert them to MRIaggr
param <- gsub(".nii", "", param)

MRIaggr.Pat1 <- constMRIaggr(ls.array, identifier = "Pat1", param = param)

#### additionnal examples
\dontrun{
## load an analyse file (example of oro.nifti::readANALYZE)
path.Pat2 <- system.file("anlz", package = "oro.nifti")
analyse.avg <- readMRI(file.path(path.Pat2, "avg152T1"), format = "analyze")
MRIaggr.Pat2 <- constMRIaggr(analyse.avg, param = "avg", identifier = "Pat2")

# display
multiplot(MRIaggr.Pat2, param = "avg",
          mfrow = c(4,6), axes = FALSE, main = "",
          mar = c(0,0,0.75,0), mar.legend = c(0,0,1,0))
		  
### load a nifti file (example of oro.nifti::readNIfTI)
path.Pat3 <- system.file("nifti", package = "oro.nifti")
nifti.ffd <- readMRI(file.path(path.Pat3, "filtered_func_data"), format = "nifti")
MRIaggr.Pat3 <- constMRIaggr(lapply(1:dim(nifti.ffd)[4], 
                                    function(x){nifti.ffd[,,,x, drop = FALSE]}),
             param=paste("ffd", 1:dim(nifti.ffd)[4], sep = "_"), identifier = "Pat3")

# display
multiplot(MRIaggr.Pat3, param = "ffd_1")

## load a dicom file (examples of oro.dicom::readDICOMFile)
path.Pat4 <- system.file("dcm", package = "oro.dicom")
dicom.Abdo <- readMRI(file.path(path.Pat4, "Abdo.dcm"), format = "dicom")
path.Pat4 <- constMRIaggr(dicom.Abdo, param = "Abdo",identifier = "Pat4")

# display
multiplot(path.Pat4, param = "Abdo")
multiplot(path.Pat4, xlim = c(100, 200), param = "Abdo")

}
}
\concept{const.}
\keyword{functions}

