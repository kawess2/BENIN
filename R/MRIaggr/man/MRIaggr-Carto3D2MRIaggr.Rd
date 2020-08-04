\name{Carto3D2MRIaggr}
\title{Carto3D to MRIaggr converter}
\alias{Carto3D2MRIaggr}
\description{
Construct a \code{\linkS4class{MRIaggr}} object by aggregating \code{\linkS4class{Carto3D}} objects.
}
\usage{
Carto3D2MRIaggr(ls.Carto3D, rm.Carto3D = FALSE, tol = 10^{-10},
         num = NULL, verbose = optionsMRIaggr("verbose"))
}
\arguments{
  \item{ls.Carto3D}{a list of \code{\linkS4class{Carto3D}} objects. REQUIRED.}
  \item{rm.Carto3D}{should the object on which the \code{ls.Carto3D} argument points be removed form the global environment ? \emph{logical}.}
  \item{tol}{numeric precision for the consistency check. \emph{positive numeric}.}
  \item{num}{the slices to extract. \emph{numeric vector} or \code{NULL}.}
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
}
\details{
ARGUMENTS: \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.
}
\value{
  a \code{\linkS4class{MRIaggr}} object.
}
\examples{
## load NIFTI files 
path.Pat1 <- system.file(file.path("nifti"), package = "MRIaggr")

Pat1.TTP.t0.nifti <- readMRI(file.path(path.Pat1, "TTP_t0"), format = "nifti")
Pat1.DWI.t0.nifti <- readMRI(file.path(path.Pat1, "DWI_t0"), format = "nifti")
Pat1.MASK_DWI.t0.nifti <- readMRI(file.path(path.Pat1, "MASK_DWI_t0"), format = "nifti")
Pat1.MASK_T2_FLAIR.t2.nifti <- readMRI(file.path(path.Pat1, "MASK_T2_FLAIR_t2"),
                                       format = "nifti")

## convert them to Carto3D
Pat1.TTP.t0.Carto3D <- constCarto3D(Pat1.TTP.t0.nifti,
             identifier = "Pat1", param = "TTP_t0", default_value = NA)
Pat1.DWI.t0.Carto3D <- constCarto3D(Pat1.DWI.t0.nifti,
             identifier = "Pat1", param = "DWI_t0", default_value = NA)
Pat1.MASK_DWI.t0.Carto3D <- constCarto3D(Pat1.MASK_DWI.t0.nifti,
             identifier = "Pat1", param = "MASK_DWI_t0", default_value = NA)
Pat1.MASK_T2_FLAIR.t2.Carto3D <- constCarto3D(Pat1.MASK_T2_FLAIR.t2.nifti,
             identifier = "Pat1", param = "MASK_T2_t2", default_value = NA)

## convert Carto3D to MRIaggr								 
MRIaggr.Pat1 <- Carto3D2MRIaggr(list(Pat1.TTP.t0.Carto3D,
                                     Pat1.DWI.t0.Carto3D,
                                     Pat1.MASK_DWI.t0.Carto3D,
                                     Pat1.MASK_T2_FLAIR.t2.Carto3D)
)
}
\concept{convert.}
\keyword{functions}
