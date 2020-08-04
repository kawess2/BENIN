\name{selectDefault_value}
\title{Extract reference values}
\alias{selectDefault_value}
\alias{selectDefault_value,Carto3D-method}
\alias{selectDefault_value,MRIaggr-method}
\usage{
\S4method{selectDefault_value}{Carto3D}(object)

\S4method{selectDefault_value}{MRIaggr}(object, param = NULL, character2numeric = FALSE)
}
\description{
Extract the reference values of the contrast parameters from a \code{\linkS4class{Carto3D}} or from a \code{\linkS4class{MRIaggr}} object.
}
\arguments{
  \item{object}{an \code{object} of class \code{\linkS4class{Carto3D}} or \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{param}{the contrast parameters for which the reference values should be returned. \emph{character vector} or \code{NULL} leading to extract reference values for all available contrast parameters.}
  \item{character2numeric}{should the default values be converted from character to numeric. \emph{logical}.}
}
\value{
  A \emph{data.frame} with one line an several columns containing the default value of each parameter.
}
\examples{
#### 1- Carto3D method ####
## load nifti files and convert them to Carto3D
path.Pat1 <- system.file("nifti", package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file.path(path.Pat1, "TTP_t0"), format = "nifti")
Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0, identifier = "Pat1", param = "TTP_t0")

## selection
selectDefault_value(Carto3D.Pat1_TTP_t0)

#### 2- MRIaggr method ####
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## select default values for all parameters
res <- selectDefault_value(MRIaggr.Pat1_red)

## select default values for T2_FLAIR_t2 and DWI_t0
res <- selectDefault_value(MRIaggr.Pat1_red, param=c("T2_FLAIR_t2", "DWI_t0"))

## select default values for T2_FLAIR_t2 and DWI_t0 and convert them in numeric
res <- selectDefault_value(MRIaggr.Pat1_red, param=c("T2_FLAIR_t2", "DWI_t0"), 
         character2numeric = TRUE)
}
\concept{select.}
\keyword{methods}

