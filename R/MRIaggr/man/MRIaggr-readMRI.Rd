\name{readMRI}
\title{Read an image file}
\alias{readMRI}
\description{
Read an image file and convert it into an array.
}
\usage{
readMRI(filename, format, na.value = 0, 
		what = "numeric", size = "NA_integer_", dimensions = NULL, 
		SPM = FALSE, reorient = FALSE, flipud = FALSE)
}
\arguments{
  \item{filename}{the file name of the image file. \emph{character}. REQUIRED.}
  \item{format}{the format of the image file. Can be \code{"raw.gz"}, \code{"analyze"},  \code{"nifti"} or \code{"dicom"}. REQUIRED.} 
  \item{na.value}{the value with which NA values are replaced. \emph{numeric} or NA.}    
  \item{what}{an object whose mode will give the mode of the vector to be read, or a character vector of length one describing the mode: one of "numeric", "double", "integer", "int", "logical", "complex", "character", "raw". Only active if \code{format} equals \code{rawb.gz}.}
  \item{size}{the number of bytes per element in the byte stream. \emph{integer}. See the documentation of the \code{writeBin} function for more details.}
  \item{dimensions}{the number of bytes per element in the byte stream. The default, NA_integer_, uses the natural size. Size changing is not supported for raw and complex vectors. Only active if \code{format} equals \code{rawb.gz}.}
  \item{SPM}{is a logical variable (default = FALSE) that forces the voxel data values to be rescaled using the funused1 ANALYZE header field. This is an undocumented convention of ANALYZE files processed using the Statistical Parametric Mapping (SPM) software. Only active if \code{format} equals \code{analyse}.}
  \item{reorient}{is a logical variable (default = TRUE) that enforces Qform/Sform transformations. Only active if \code{format} equals \code{nifti}.}
  \item{flipud}{is a logical variable for vertical flipping of the image (default is TRUE). Only active if \code{format} equals \code{dicom}.}
  }
\details{

This function requires to have installed the \emph{oro.nifti} package to work if argument format is set to \code{"analyze"} or \code{"nifti"}.
It requires to have installed the \emph{oro.dicom} package to work if argument format is set to dicom.

ARGUMENTS :
\code{filename} argument corresponds to:
\itemize{
\item the \code{con} argument of the \code{base::readBin} function. 
\item the \code{fname} argument of the \code{oro.nifti::readANALYZE} function. It should be a pathname to .img or .hdr files without the suffix.
\item the \code{fname} argument of the \code{oro.nifti::readNIfTI} function. 
\item the \code{fname} argument of the \code{oro.dicom::readDICOMFile} function. It should be the file name with suffix. 
}
The \code{what} and \code{dimensions} correspond to the \code{what} and \code{size}  argument of the \code{base::readBin} function. \cr
The \code{SPM} corresponds to the \code{SPM} argument of the \code{oro.nifti::readANALYZE} function. \cr
The \code{reorient} corresponds to the \code{reorient} argument of the \code{oro.nifti::readNIfTI} function. \cr
The \code{flipud} corresponds to the \code{flipud} argument of the \code{oro.dicom::readDICOMFile} function. \cr

FUNCTION : \cr
This function is a copy of the \code{readMRI} function of the mritc package. It calls
\code{readANALYZE} of the oro.nifti package to read analyze files, 
\code{readNIfTI} of the oro.nifti package to read NIFTI files, 
\code{readDICOMFile} of the oro.dicom package to read dicom files 
and \code{readBin} of the base package to read rawb.gz.
}
\value{
  A \code{nifti} object if \code{format="nifit"}, an \code{anlz} object if \code{format="analyse"} and else an \code{array}.
}
\examples{ 
## load a NIFTI file
path.Pat1 <- system.file(file.path("nifti"), package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file.path(path.Pat1, "TTP_t0"), format = "nifti")
dim(nifti.Pat1_TTP_t0)

\dontrun{
## load an analyse file (example of oro.nifti::readANALYZE)
path.Pat2 <- system.file("anlz", package = "oro.nifti")
analyse.avg <- readMRI(file.path(path.Pat2, "avg152T1"), format = "analyze")
graphics::image(analyse.avg[,,45,1])

### load a NIFTI file (example of oro.nifti::readNIfTI)
path.Pat3 <- system.file("nifti", package = "oro.nifti")
nifti.ffd <- readMRI(file.path(path.Pat3, "filtered_func_data"), format = "nifti")
graphics::image(nifti.ffd[,,10,32])

## load a dicom file (examples of oro.dicom::readDICOMFile)
path.Pat4 <- system.file("dcm", package = "oro.dicom")
dicom.Abdo <- readMRI(file.path(path.Pat4,"Abdo.dcm"), format = "dicom")
graphics::image(dicom.Abdo[[2]], col = grey(0:64/64), axes = FALSE, xlab = "", ylab = "",
                main = "Abdo.dcm")
}
}
\concept{convert.}
\keyword{functions}
