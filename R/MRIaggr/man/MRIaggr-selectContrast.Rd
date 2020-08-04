\name{selectContrast}
\title{Extract contrast parameters}
\alias{selectContrast}
\alias{selectContrast,Carto3D-method}
\alias{selectContrast,MRIaggr-method}
\description{
Extract the contrast parameters from a \code{\linkS4class{Carto3D}} or a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{selectContrast}{Carto3D}(object, num = NULL, na.rm = FALSE, coords = FALSE,
         format = "data.frame")
			
\S4method{selectContrast}{MRIaggr}(object, param = NULL, num = NULL, format = "data.frame",
         slice_var = "k", coords = FALSE, hemisphere = "both", 
		 norm_mu = FALSE, norm_sigma = FALSE, na.rm = FALSE, subset = NULL)
}
\arguments{
  \item{object}{an \code{object} of class \code{\linkS4class{Carto3D}} or \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{param}{the contrast parameters to extract. \emph{character vector} or \code{NULL}.}
  \item{num}{the slices to extract. \emph{numeric vector} or \code{NULL}.}
  \item{format}{the format of the output. Can be \code{"matrix"}, \code{"data.frame"} or \code{"vector"}.}
  \item{slice_var}{the type of slice to extract. \code{"i"} for sagittal, \code{"j"} for coronal and \code{"k"} for transverse. \emph{character}.}  
  \item{coords}{the coordinates that sould be extracted. \emph{logical} or any of \code{"i"} \code{"j"} \code{"k"}.}
  \item{hemisphere}{the hemisphere to extract. \emph{character}. See the details section.}
  \item{norm_mu}{the type of centering to apply on the parameter values. \emph{character}. See the details section.}
  \item{norm_sigma}{the type of scaling to apply on the parameter values. \emph{character}. See the details section.}
  \item{na.rm}{should observations with missing values be removed ? \emph{logical}.}
  \item{subset}{the subset of observations to extract. \emph{positive integer vector} or \code{NULL} leading to use all observations}
}
\details{
ARGUMENTS: \cr
Information about the \code{param} argument can be found in the details section of \code{\link{initParameter}}.

Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

Possible values for the \code{hemisphere} argument are:
\itemize{
\item \code{"both"} : select all the observations.
\item \code{"left"} : select the observations from the left hemisphere.
\item \code{"right"} : select the observations from the right hemisphere.
\item \code{"lesion"} : select the observations belonging to the hemisphere(s) that contain(s) the lesion (if any). 
\item \code{"contralateral"} : select the observations belonging to the hemisphere(s) that do not contain(s) the lesion (if any).
}
To select observations from a given hemisphere (all values except \code{"both"}), the parameter \code{hemisphere} must have been allocated to the object using, for instance, \code{\link{calcHemisphere}}.
In addition for \code{"lesion"} and \code{"contralateral"} values,  the slot \code{@hemispheres} has to be filled using, for instance, \code{\link{calcHemisphere}}.

Possible values for the centering argument (\code{norm_mu}) and the scaling argument (\code{norm_sigma}) are:
\itemize{
\item \code{"FALSE"} : no normalization
\item \code{"global"} : the centering or scaling value is computed using all the observations.
\item \code{"global_1slice"} : the centering or scaling value is computed using all the observations that belong to the slice of the observation to normalize.
\item \code{"global_3slices"} : the centering or scaling value is computed using all the observations that belong to the slice of the observation to normalize, the slice above (if any) and the slice below (if any).
\item \code{"contralateral"} : the centering or scaling value is computed using the observations from the contralateral hemisphere.
\item \code{"contralateral_1slice"} :  the centering or scaling value is computed using the observations from the contralateral hemisphere that belong to the slice of the observation to normalize.
\item \code{"contralateral_3slices"} :  the centering or scaling value is computed using the observations from the contralateral hemisphere that belong to the slice of the observation to normalize, the slice above (if any) and the slice below (if any).
\item \code{"default_value"} : the default value of the parameter stored in the slot \code{@default_value} is used for the centering (for \code{norm_mu} only).
}

If \code{coords} is set to \code{TRUE} the dataset containing the contrast parameters values will also contains all the coordinates. If \code{coords} is set to \code{FALSE}, it will not contain any coordinates. 

Argument \code{subset} can be a \emph{character} value that refers to a logical parameter in the \code{object} defining the subset of observation to extract.

FUNCTION: \cr
Each of the \code{num}, \code{hemisphere} and \code{subset} argument define a subset of the total set of observations. 
It is the intersection of all these three subsets that is extracted.

When a normalisation is requested to center (resp. scale) the data, the normalisation value is extracted for each parameter in the element of the slot normalization that match the argument \code{norm_mu} (resp. \code{norm_sigma}).
The parameters values are first centered by substraction with the value returned by \code{norm_mu}. 
Then they are scaled by division with the value returned by \code{norm_sigma}.

}
\value{
   A \code{"vector"}, \emph{matrix} or a \emph{data.frame}. In the latter two cases, each row corresponds to a voxel and each column to a coordinate.
}
\seealso{ 
 \code{\link{calcContralateral}}, \code{\link{calcRegionalContrast}}, \code{\link{calcFilter}} and \code{\link{calcTissueType}} to retreat and allocate the modified contrast parameters. \cr
 \code{\link{allocContrast<-}} to allocate new contrast parameters. \cr 
 \code{\link{calcNormalization}} to compute and allocate the normalisation values. \cr
 \code{\link{allocNormalization<-}} to allocate the normalization values when obtained from an external source. \cr
 \code{\link{calcHemisphere}} and \code{\link{calcContralateral}} to compute and allocate the hemispheres. \cr
 \code{\link{allocHemisphere<-}} and \code{\link{allocContrast<-}} to allocate hemispheres obtained from an external source.
}
\examples{
#### 1- Carto3D method ####
## load nifti files and convert them to Carto3D
path.Pat1 <- system.file("nifti", package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file.path(path.Pat1, "TTP_t0"), format = "nifti")
Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0, identifier = "Pat1", param = "TTP_t0")

## select all observations
carto1 <- selectContrast(Carto3D.Pat1_TTP_t0)
dim(carto1)

## select observations from slices 1 to 3 and return the result into a data.frame
carto2 <- selectContrast(Carto3D.Pat1_TTP_t0, num = 1:3, coords = FALSE, format = "data.frame")
dim(carto2)

## select observations from slices 1 to 3 and return the result into a vector
carto3 <- selectContrast(Carto3D.Pat1_TTP_t0, num = 1:3, coords = FALSE)
length(carto3)

#### 2- MRIaggr method ####
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## select all parameters and all observations
carto <- selectContrast(MRIaggr.Pat1_red)
dim(carto)
head(carto)

## select a subset of parameters
carto <- selectContrast(MRIaggr.Pat1_red, param = c("DWI_t0","T2_FLAIR_t2"))
dim(carto)
head(carto)

## select a subset of parameters on slices 1 to 3
carto <- selectContrast(MRIaggr.Pat1_red, num=1:3, param=c("DWI_t0","T2_FLAIR_t2"))
dim(carto)
head(carto)

## select a subset of parameters on slices 1 to 3 and normalized the center 
## the values using the contralateral
carto <- selectContrast(MRIaggr.Pat1_red, num=1:3, param=c("DWI_t0","T2_FLAIR_t2"),
                        norm_mu = "contralateral")
dim(carto)
head(carto)

## select only observations which are lesioned at admission (i.e. MASK_DWI_t0 = TRUE)
carto <- selectContrast(MRIaggr.Pat1_red, subset = "MASK_DWI_t0",
                        param = c("DWI_t0","T2_FLAIR_t2","MASK_DWI_t0"))
dim(carto)
head(carto)

## select only observations which are lesioned at admission (i.e. MASK_DWI_t0 = TRUE) 
## with coordinates
carto <- selectContrast(MRIaggr.Pat1_red, subset = "MASK_DWI_t0",
                        param = c("DWI_t0","T2_FLAIR_t2","MASK_DWI_t0"), coords = TRUE)
dim(carto)
head(carto)

## select only observations for which i = 55
carto <- selectContrast(MRIaggr.Pat1_red, slice_var = "i", 
                        num = 55, coords = TRUE)
dim(carto)
head(carto)
}
\concept{select.}
\keyword{methods}

