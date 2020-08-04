\name{array2df}
\title{Array to data.frame converter}
\alias{array2df}
\description{
Convert observations stored in the array format into the data.frame format.
}
\usage{
array2df(array, coords = NULL, name_newparam = "res",
         names_coords = letters[9:(8+ncol(coords))], na.rm = TRUE)
}
\arguments{
  \item{array}{the array that should be converted into a data.frame. \emph{array} or \emph{matrix}. REQUIRED.}
  \item{coords}{the spatial coordinates of the observations contained in \code{array}. \emph{matrix} or \code{NULL}.}
  \item{name_newparam}{the name of the contrast parameter to which cooresponds \code{array}. \emph{character}.}
  \item{names_coords}{the name of the coordinates. \emph{character vector}.}
  \item{na.rm}{should observations with missing values be removed ? \emph{logical}.}
}
\details{
ARGUMENTS : \cr
If \code{coords} is set to \code{NULL}, the coordinates will be defined by the position of the observations in \code{array}. \cr
If \code{na.rm} is set to \code{TRUE}, the \code{coord} argument will only contains the coordinates of the non-NA observations of \code{array}.
}
\value{
  A \code{dataframe} with in columns the coordinates and the parameter values, and in rows the observations.
}
\examples{ 
## load a NIFTI file (array format)
path.Pat1 <- system.file(file.path("nifti"), package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file.path(path.Pat1, "TTP_t0"), format = "nifti")
dim(nifti.Pat1_TTP_t0)

## conversion to data frame format
res128 <- array2df(array = nifti.Pat1_TTP_t0, name_newparam = "TTP_t0")
dim(res128)
head(res128)

## conversion to data frame format with specific coordinates
res256 <- array2df(array = nifti.Pat1_TTP_t0, name_newparam = "TTP_t0",
                   coords = expand.grid(129:206, 129:228, 1:3, 1))
dim(res256)
head(res256)

}
\concept{convert.}
\keyword{functions}
