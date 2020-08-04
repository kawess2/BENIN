\name{allocContrast}
\title{Allocate new contrast parameters}
\alias{allocContrast<-}
\alias{allocContrast<-,MRIaggr-method}
\description{
Allocate one or several contrast parameters to a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{allocContrast}{MRIaggr}(object, param = NULL, default_value = NULL,
         overwrite = FALSE, verbose = optionsMRIaggr("verbose")) <- value
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{value}{the value of each contrast parameter (in columns) at each voxel (in rows). \emph{data.frame}. REQUIRED.}
  \item{param}{the names of the contrast parameters. \emph{character vector} or \code{NULL} leading to use the name of the \code{value} argument.}
  \item{default_value}{the reference values of the contrast parameters (e.g background values). A one row \emph{data.frame} where the column length must match the length of the \code{param} argument.}
  \item{overwrite}{if a contrast parameters with the same names are already stored in \code{object@data}, can they be overwritten ? \emph{logical}.}  
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
}
\details{
FUNCTION: \cr
If the \code{param} argument is not specified then the names of the \code{value} argument will be used to define the parameter names.

If the \code{default_value} argument is  \code{NULL} then default values \code{"undefined"} are attributed to each parameter. 

If a parameter named \code{"mask"} is intended to be allocated it must be done alone and it must be of type logical.

Parameter names \code{"index"}, \code{"i"}, \code{"j"} and \code{"k"} are reserved and cannot be modified.
Nevertheless if \code{value} contains parameters \code{"i"}, \code{"j"} and \code{"k"}, the correspondance between these coordinates and the object coordinates is tested.
}
\seealso{
  \code{\link{calcContralateral}}, \code{\link{calcRegionalContrast}}, \code{\link{calcFilter}} and \code{\link{calcTissueType}} to compute, modify and allocate cartography. \cr
  \code{\link{selectContrast}} to select contrast parameters in the \code{MRIaggr} object. 
}
\value{
  None.
}
\examples{
## load NIFTI files and convert them to MRIaggr
path.Pat1 <- system.file(file.path("nifti"), package = "MRIaggr")
ls.array <- list()
ls.array[[1]] <- readMRI(file.path(path.Pat1, "DWI_t0"),format = "nifti")
ls.array[[2]] <- readMRI(file.path(path.Pat1, "MASK_DWI_t0"),format = "nifti")
MRIaggr.Pat1 <- constMRIaggr(ls.array, identifier = "Pat1", param = c("DWI_t0", "MASK_DWI_t0"))

## allocate a new contrast parameters
allocContrast(MRIaggr.Pat1, param = "noise", overwrite = TRUE) <- rnorm(selectN(MRIaggr.Pat1))

## perform operations on a contrast parameters and store the results
myCarto <- selectContrast(MRIaggr.Pat1 , param = "DWI_t0")
myCarto <- myCarto * 2 + 1
allocContrast(MRIaggr.Pat1, param = "myCarto", overwrite = TRUE) <- myCarto

## import a contrast parameters in an already existing MRIaggr object
nifti.MTT_t0 <- readMRI(file.path(path.Pat1, "MTT_t0"), format = "nifti")
df.MTT_t0 <- array2df(nifti.MTT_t0, name_newparam = "MTT_t0")$MTT_t0
allocContrast(MRIaggr.Pat1, param = "MTT_t0", overwrite = TRUE) <- df.MTT_t0

## some calc methods automatically save results in the @data slot
calcFilter(MRIaggr.Pat1, param = "MTT_t0", filter = "2D_G3",
           update.object = TRUE, overwrite = TRUE)
res <- selectContrast(MRIaggr.Pat1, param = "MTT_t0_2D_G3")

}
\concept{alloc.}
\keyword{methods}

