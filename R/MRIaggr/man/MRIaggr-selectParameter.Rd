\name{selectParameter}
\title{Extract parameters}
\alias{selectParameter}
\alias{selectParameter,Carto3D-method}
\alias{selectParameter,MRIaggr-method}
\usage{
\S4method{selectParameter}{Carto3D}(object)

\S4method{selectParameter}{MRIaggr}(object, type = "contrast", mask = TRUE)
}
\description{
Extract parameters from a \code{\linkS4class{Carto3D}} or from a \code{\linkS4class{MRIaggr}} object.
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{type}{the type of parameter to select. \emph{character}. See the details section.}
  \item{mask}{should the mask be considered as an available contrast parameter ? \emph{logical}.}
}
\details{
ARGUMENTS: \cr
Possible values for \code{type} are:
\itemize{
\item \code{"clinic"} : return the name of the clinical parameters.
\item \code{"contrast"} : return the name of the contrast parameters.
\item \code{"ls_descStats"} : return the name of the elements in \code{@ls_descStats}.
}
}
\seealso{
 \code{\link{allocContrast<-}} to allocate new parameters. 
}
\examples{
#### 1- Carto3D method ####
## load NIFTI files and convert them to Carto3D
path.Pat1 <- system.file("nifti", package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file.path(path.Pat1, "TTP_t0"), format = "nifti")
Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0, identifier = "Pat1", param = "TTP_t0")

## selection
selectParameter(Carto3D.Pat1_TTP_t0)

#### 2- MRIaggr method ####
## load a MRIaggr object
data("MRIaggr.Pat1_red", package="MRIaggr")

## extract all imaging parameters
res <- selectParameter(MRIaggr.Pat1_red)

## extract the names of the parameters in the slot @ls_descStats
res <- selectParameter(MRIaggr.Pat1_red, type = "ls_descStats")

}
\concept{alloc.}
\keyword{methods}
