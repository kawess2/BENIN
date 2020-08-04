\name{multiplot}
\title{Slice by slice display}
\alias{multiplot}
\alias{multiplot,data.frame-method}
\alias{multiplot,Carto3D-method}
\alias{multiplot,MRIaggr-method}
\description{
Make a slice by slice display of a \emph{data.frame}, a \code{\linkS4class{Carto3D}} or a \code{\linkS4class{MRIaggr}} object.
}
\usage{
\S4method{multiplot}{data.frame}(object, contrast = NULL, num = NULL, 
         index1 = NULL, index2 = NULL, index3 = NULL,
         col = NULL, pch = NULL ,xlim = NULL, ylim = NULL, filename = "multiplot",...)
		 
\S4method{multiplot}{Carto3D}(object, num = NULL, 
         col=NULL, pch = NULL, xlim = NULL, ylim = NULL,
         filename = "multiplot",...)

\S4method{multiplot}{MRIaggr}(object, param, num = NULL,
         index1 = NULL, index2 = NULL, index3 = NULL,
         midplane = FALSE, col = NULL, pch = NULL, xlim = NULL, ylim = NULL,
		 filename = "multiplot", ...)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}, or a \code{\linkS4class{Carto3D}}, or a 2 or 3 column \code{data.frame} containing the coordinates of the observations in columns. REQUIRED.}
  \item{param}{the contrast parameter to display. \emph{character}. REQUIRED.}
  \item{contrast}{the intensities to display. \emph{numerical vector} or \code{NULL} leading to use the same color for all observations.}
  \item{num}{the slices to display. \emph{numeric vector} or \code{NULL}.}
  \item{index1,index2,index3}{the coordinates of additionnal points to display. \emph{data.frame} or \emph{list} or \code{NULL}.}
  \item{midplane}{should the mid-sagittal plan be displayed ? \emph{logical}.}
 
  \item{col}{the color with which the observations will be displayed. \emph{character vector} or \code{NULL} leading to determine the colors using the \code{palette} and \code{breaks} arguments.}
  \item{pch}{the symbol with which the observations will be displayed. \emph{positive integer} or \code{NULL} leading to use the \code{image} function instead of \code{plot}.}

  \item{xlim}{the x limits of the plot. \emph{numeric vector of size 2} or \code{NULL} leading to automatic setting of the x limits.}
  \item{ylim}{the y limits of the plot. \emph{numeric vector of size 2} or \code{NULL} leading to automatic setting of the y limits.}
  \item{filename}{the name of the file used to export the plot. \emph{character}.}
  \item{...}{additional arguments to be passed to \code{\link{optionsMRIaggr}} for specifying the graphical parameters.}
}
\details{
ARGUMENTS: \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.
Information about the \code{index1}, \code{index2} and \code{index3} arguments can be found in the details section of \code{\link{initIndex}} (argument \code{index}).

Information about the \code{...} can be found in the documentation of the \code{\link{optionsMRIaggr}} function.

}
\value{
  Invisible, return a list containing :
  \tabular{ll}{
  \code{breaks.plot}         \tab the breaks used to display the images. \cr
  \code{palette.plot}        \tab the palette used to display the images. \cr
  \code{breaks.legend}       \tab the breaks used to display the legend. \emph{character}. \cr
  \code{palette.legend}      \tab the palette used to display the legend. \emph{character}. \cr
  \code{quantiles.legend}    \tab the quantiles used to display the legend. \emph{character}. \cr
  }
}

\seealso{ 
 \code{\link{plotLesion3D}} for a 3D plot of the lesion. \cr
 \code{slices3d} of the \emph{misc3d} package for a more interactive 3D plot.
}
\examples{
#### 1- data.frame ####
## simulate 
n <- 10
Y <- rnorm(n^2)

## display
multiplot(object = data.frame(expand.grid(1:n,1:n),1),
             contrast = Y, window = FALSE)

## load a MRIaggr object
data(MRIaggr.Pat1_red, package = "MRIaggr")

## select data
data <- selectContrast(MRIaggr.Pat1_red, 
         param = c("DWI_t0","TTP_t0","MTT_t0","MASK_T2_FLAIR_t2"),
         hemisphere = "lesion", coords = TRUE)
					
## fit model
glm.1 <- glm(MASK_T2_FLAIR_t2 ~ DWI_t0 + TTP_t0 + MTT_t0, data = data, 
         family = binomial(link = "logit"))

## display fitted values
multiplot(object = data[,c("i","j","k")],
             contrast = predict(glm.1, type = "response"), window = FALSE)

## display residuals
multiplot(object = data[,c("i","j","k")], num = 3,
             contrast = predict(glm.1, type = "response"), window = FALSE,
             index1 = list(coords = data[data$MASK_T2_FLAIR_t2,c("i","j","k")], outline = TRUE)
             )
			  
#### 2- carto3D ####
## load NIFTI files and convert them to carto3D
path.Pat1 <- system.file("nifti", package = "MRIaggr")
nifti.Pat1_TTP_t0 <- readMRI(file.path(path.Pat1, "TTP_t0"), format = "nifti")
Carto3D.Pat1_TTP_t0 <- constCarto3D(nifti.Pat1_TTP_t0, identifier = "Pat1", param = "TTP_t0")

## display
multiplot(Carto3D.Pat1_TTP_t0)
\dontrun{
multiplot(Carto3D.Pat1_TTP_t0, num = 1:2)
multiplot(Carto3D.Pat1_TTP_t0, num = 1:2, axes = FALSE)
multiplot(Carto3D.Pat1_TTP_t0, num = 1:2, axes = FALSE, legend = FALSE)
multiplot(Carto3D.Pat1_TTP_t0, num = 1:2, axes = FALSE, legend = FALSE, 
         main = "", num.main = FALSE)
multiplot(Carto3D.Pat1_TTP_t0, num = 1:2, axes = FALSE, main = "", num.main = FALSE,
			 palette = "gray.colors", breaks = seq(0,100))
}

#### 3- MRIaggr ####
## load a MRIaggr object
data(MRIaggr.Pat1_red, package = "MRIaggr")

# display 3 slices 
multiplot(MRIaggr.Pat1_red, param = "DWI_t0",              
             num = 1:3)

\dontrun{
# display 3 slices with no axes and white background
multiplot(MRIaggr.Pat1_red, param = "DWI_t0",
             num = 1:3, axes = FALSE, bg = "white")

# remove the legend
multiplot(MRIaggr.Pat1_red, param = "DWI_t0",              
             num = 1:3, legend = FALSE)
}

## display an set of points
# using a binary parameter stored in the object
multiplot(MRIaggr.Pat1_red, param = "DWI_t0",              
             num = 1:3, index1 = list(coords = "MASK_DWI_t0")
)

\dontrun{
# customize the display of the points
multiplot(MRIaggr.Pat1_red, param = "DWI_t0",              
             num = 1:3, index1 = list(coords = "MASK_DWI_t0", col = "pink", pch = 14)
)

# display only the edges of the set
multiplot(MRIaggr.Pat1_red, param = "DWI_t0", num = 3, legend = FALSE,  
             index1 = list(coords = "MASK_DWI_t0", outline = TRUE)
)
             

# specify the index of points using coordinates
coordsIndex <- data.frame(i = c(40,60), j = c(80,100), k = c(3,3))

multiplot(MRIaggr.Pat1_red, param = "DWI_t0", num = 3, legend = FALSE,       
             index2 = list(coords = coordsIndex, col = "black", pch = 15, cex = 4)
             )
			 
# various possibilities for the display
multiplot(MRIaggr.Pat1_red, num = 1:3, param = "DWI_t0",
             legend = FALSE, window = FALSE)
multiplot(MRIaggr.Pat1_red, num = 1:3, param = "DWI_t0",
             legend = TRUE, window = FALSE)
multiplot(MRIaggr.Pat1_red, num = 1:3, param = "DWI_t0",
             legend = NULL, window = FALSE)
multiplot(MRIaggr.Pat1_red, num = 1:3, param = "DWI_t0",
             legend = "only", window = FALSE)
			 }
}
\concept{plot.}
\keyword{methods}

