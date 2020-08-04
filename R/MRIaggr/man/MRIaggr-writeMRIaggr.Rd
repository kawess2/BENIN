\name{writeMRIaggr}
\title{Write an image file}
\alias{writeMRIaggr}
\alias{writeMRIaggr,MRIaggr-method}
\description{
Write an image file from a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{writeMRIaggr}{MRIaggr}(object, param, num = NULL, 
         norm_mu = FALSE, norm_sigma = FALSE, range.coords = NULL, 
		 default_value = NA, filename, format, gzipped = TRUE,
         verbose = optionsMRIaggr("verbose"), size = "NA_integer_")
}
\arguments{
	\item{object}{an \code{object} of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
	\item{param}{the contrast parameters to extract. \emph{character vector} or \code{NULL}. REQUIRED.}
	\item{num}{the slices to extract. \emph{numeric vector} or \code{NULL}.}
	\item{norm_mu}{the type of centering to apply on the parameter values. \emph{character}. See the details section.}
	\item{norm_sigma}{the type of scaling to apply on the parameter values. \emph{character}. See the details section.}
	\item{default_value}{the element used to fill the missing observations. \emph{numeric}.}
    \item{range.coords}{the maximum coordinate in each dimension to be considered. \emph{numeric vector} of length 3.}
	\item{filename}{the name of the image file. \emph{character}. REQUIRED.}
	\item{format}{the format of the image file. Can be \code{"raw.gz"}, \code{"analyze"},  \code{"nifti"} or \code{"dicom"}. REQUIRED.} 
	\item{size}{the number of bytes per element in the byte stream. \emph{integer}. See the documentation of the \code{writeBin} function for more details.}
	\item{gzipped}{is a character string that enables exportation of compressed (.gz) files.}
	\item{verbose}{is a logical variable that allows text-based feedback during execution of the function.} 
 }
\details{
This function requires to have installed the \emph{oro.nifti} package to work if argument format is set to \code{"analyze"} or \code{"nifti"}.
}
\value{
  None.
}
\examples{ 
data("MRIaggr.Pat1_red", package = "MRIaggr")

\dontrun{
writeMRIaggr(MRIaggr.Pat1_red, param = "DWI_t0", filename = "Pat1_DWIred", format = "nifti")
DWIred <- readMRI("Pat1_DWIred", format = "nifti")
image(DWIred[,,1])
}

}
\concept{convert.}
\keyword{methods}
