\name{calcFilter}
\title{Image filtration}
\alias{calcFilter}
\alias{calcFilter,array-method}
\alias{calcFilter,MRIaggr-method}
\description{
Apply a filter to an image.
}
\usage{
\S4method{calcFilter}{array}(object, filter, norm.filter = TRUE, bilateral = FALSE, na.rm = FALSE)

\S4method{calcFilter}{MRIaggr}(object, param, filter, norm.filter = TRUE, bilateral = FALSE,
     na.rm = FALSE, name_newparam = NULL,
	 verbose = optionsMRIaggr("verbose"), update.object = FALSE, overwrite = FALSE)
}
\arguments{
  \item{object}{an \code{array} or an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{param}{the contrast parameter to be filtered. \emph{character vector}. REQUIRED.}
  \item{filter}{the filter to use. Can be a \emph{matrix} or an \emph{array}, or a name indicating which filter should be used. REQUIRED.}
  \item{norm.filter}{should the filtered correspond to a weighted mean over site ? (or a weighted sum). \emph{logical}.}
  \item{bilateral}{should the influence of each neighbor be ponderated by the difference in signal with the considered observation ? \emph{logical}.}
  \item{na.rm}{should observations with missing values in their neighbourhood be set to NA ? Otherwise the ponderation is adjusted. \emph{logical}.}
  \item{name_newparam}{the name of the new parameters. \emph{character vector}.}
  
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
  \item{update.object}{should the resulting filtered parameters be stored in \code{object} ? \emph{logical}.}
  \item{overwrite}{if contrast parameters with the same names are already stored in \code{object@data}, can they be overwritten ? \emph{logical}.}  
}
\details{
ARGUMENTS: \cr
Several types of pre-stored filters are availables and can be called by their name: 
\itemize{
\item filter for smoothing purpose (e.g. "3D_G5") : see the details section of \code{\link{initFilter}}.
\item filter for neighbourhood definition purpose (e.g. "3D_N10") : see the details section of \code{\link{initNeighborhood}}.
}
\code{norm.filter} should be set to \code{TRUE} for gaussian, median or sobel filters to correct edge effects. 
This lead to weight the neighboring value in order to offset the incomplete neighbourhood.
}
\seealso{
\code{\link{selectContrast}} to select the filtered parameter(s).
\code{\link{initFilter}} or code{\link{initNeighborhood}} to select pre-stored filters.
}
\value{
  An \emph{list} containing :
  \itemize{
  \item\code{[[res]]} : a \emph{data.frame} containing the coordinates and the filtered parameters.
  \item\code{[[filter]]} : the name of the filter that has been used. \emph{character}.
  }
}
\examples{
#### 1- array method ####
## load a NIFTI file
path.Pat1 <- system.file(file.path("nifti"), package = "MRIaggr")
nifti.Pat1_DWI_t0 <- readMRI(file.path(path.Pat1, "DWI_t0"), format = "nifti")

## before filtering
graphics::image(nifti.Pat1_DWI_t0[,,1,1])

## after median filtering
niftiF.Pat1_DWI_t0 <- calcFilter(nifti.Pat1_DWI_t0[,,,1], filter = "2D_M3")$res
graphics::image(niftiF.Pat1_DWI_t0[,,1])

#### 2- MRIaggr method ####
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## compute and allocate filtered parameter to the MRIaggr object
# gaussian filter
calcFilter(MRIaggr.Pat1_red, param = c("T2_FLAIR_t2","DWI_t0","TTP_t0"),
           filter = "2D_G3", bilateral = FALSE, na.rm = FALSE,
		   update.object = TRUE, overwrite = TRUE)
selectParameter(MRIaggr.Pat1_red)

# median filter
calcFilter(MRIaggr.Pat1_red, param = c("T2_FLAIR_t2","DWI_t0","TTP_t0"),
           filter = "2D_M3", na.rm = FALSE, update.object = TRUE, overwrite = TRUE)

## display
par(mfrow = c(2,2), mar = c(2,2,2,2))
multiplot(MRIaggr.Pat1_red, param = "T2_FLAIR_t2",
          num.main = FALSE, main = "raw",
          num = 1, window = NULL, breaks = c(-100,seq(0, 450),601))
multiplot(MRIaggr.Pat1_red, param = "T2_FLAIR_t2_2D_G3",
          num.main = FALSE, main = "2D_G3",
          num = 1, legend = FALSE, window = NULL, breaks = c(-100,seq(0,450),601))
multiplot(MRIaggr.Pat1_red, param = "T2_FLAIR_t2_2D_M3",
          num.main = FALSE, main = "2D_M3",
          num = 1, legend = FALSE, window = NULL, breaks = c(-100,seq(0, 450),601))

## see the results of the different filters
# G : Gaussian filter
resG <- calcFilter(MRIaggr.Pat1_red, param = "T2_FLAIR_t2",
                   filter = "2D_G3")

# M : median filter
resM <- calcFilter(MRIaggr.Pat1_red, param = "T2_FLAIR_t2",
                   filter = "2D_M3")

# S : Sobel filter
resS <- calcFilter(MRIaggr.Pat1_red, param = "T2_FLAIR_t2",
                  filter = "2D_Sx")

# I
resI.T <- calcFilter(MRIaggr.Pat1_red, param = "MASK_T2_FLAIR_t2",
                  filter = "2D_I3", norm.filter = TRUE)
resI.F <- calcFilter(MRIaggr.Pat1_red, param = "MASK_T2_FLAIR_t2",
                   filter = "2D_I3", norm.filter = FALSE)

# N
resN.T <- calcFilter(MRIaggr.Pat1_red, param = "MASK_T2_FLAIR_t2",
                     filter = "3D_N10", norm.filter = TRUE)
resN.F <- calcFilter(MRIaggr.Pat1_red, param = "MASK_T2_FLAIR_t2",
                     filter = "3D_N10", norm.filter = FALSE)

## display
par(mfrow = c(2,2), mar = rep(2,4), mgp = c(2,0.75,0))
breaks <- seq(-50,500,1)
multiplot(MRIaggr.Pat1_red, param = "T2_FLAIR_t2", num = 3,
             breaks = breaks, window = NULL, legend = FALSE,
             main = "no filtering", num.main = FALSE, main.legend = "")
multiplot(resS$res[,c("i","j","k")], contrast = resS$res[,4],
             num = 1, window = NULL, legend = FALSE,
             palette = "cm.colors", breaks = seq(-100, 100),
             main = "sobelX filtering", num.main = FALSE)
multiplot(resG$res[,c("i","j","k")], contrast = resG$res[,4],
             num = 3, window = NULL, legend = FALSE, breaks = breaks,
             main = "gaussian filtering", num.main = FALSE)
multiplot(resM$res[,c("i","j","k")], contrast = resM$res[,4],
             num = 3, window = NULL, legend = FALSE, breaks = breaks,
             main = "median filtering", num.main = FALSE)

layout(matrix(1:6, nrow = 3, ncol = 2, byrow = TRUE), widths = c(2,1), heights = rep(1,3))
par(mar = rep(2,4), mgp = c(2,0.75,0))
multiplot(MRIaggr.Pat1_red, param = "MASK_T2_FLAIR_t2", num = 3,
             window = NULL, legend = TRUE, main = "raw", num.main = FALSE, main.legend = "")
multiplot(resI.T$res[,c("i","j","k")], contrast = resI.T$res[,4],
             num = 3, window = NULL, legend = TRUE,
             main = "Influence filtering - norm = TRUE", num.main = FALSE)
multiplot(resI.F$res[,c("i","j","k")], contrast = resI.F$res[,4],
             num = 3, window = NULL, legend = TRUE,
             main = "Influence filtering - norm = FALSE", num.main = FALSE)


layout(matrix(1:6, nrow = 3, ncol = 2, byrow = TRUE), widths = c(2,1), heights = rep(1,3))
par(mar = rep(2,4), mgp = c(2,0.75,0))
multiplot(MRIaggr.Pat1_red, param = "MASK_T2_FLAIR_t2", num = 3,
             window = NULL, legend = TRUE, main = "raw", num.main = FALSE, main.legend = "")
multiplot(resN.T$res[,c("i","j","k")], contrast = resN.T$res[,4],
             num = 3, window = NULL, legend = TRUE,
             main = "Neighborhood3D filtering - norm = TRUE", num.main = FALSE)
multiplot(resN.F$res[,c("i","j","k")], contrast = resN.F$res[,4],
             num = 3, window = NULL, legend = TRUE,
             main = "Neighborhood3D filtering - norm = FALSE", num.main = FALSE)

}
\concept{calc.}
\keyword{methods}

