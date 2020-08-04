\name{df2array}
\title{data.frame to array converter}
\alias{df2array}
\usage{
df2array(contrast, coords, format = "any", default_value = NA, range.coords = NULL)
}
\description{
Convert observations stored in the data.frame format into the array format.
}
\arguments{
  \item{contrast}{the dataset containing the observations in rows and the contrast parameters in columns. \emph{vector} or \emph{data.frame}. REQUIRED.}
  \item{coords}{the spatial coordinates of the observations. \emph{matrix} with a number of rows equal to the number of rows of \code{data}. REQUIRED.}
  \item{format}{the format of the output. Can be \code{"any"},\code{"matrix"},\code{"data.frame"} or \code{"list"}.}
  \item{default_value}{the element used to fill the missing observations. \emph{numeric}.}
  \item{range.coords}{the maximum coordinate in each dimension to be considered. \emph{numeric vector} with length equal to the number of columns of \code{coords}.}
}
\details{
FUNCTION: \cr
If \code{contrast} contains several parameters, they are treated one at a time and the result is returned in the form of a list.
If \code{range.coords} is \code{NULL} then the maxima coordinates are those of the \code{coords} argument.
If only one parameter is specified and the format is set to \code{"any"} then a \emph{vector} is returned.
}
\value{
a list containing :
\itemize{
\item\code{[[contrast]]} :  \emph{list} containing the new contrast in the new format.
\item\code{[[coords]]} :  a \emph{data.frame} containing the coordinates of each observation.
\item\code{[[unique_coords.group]]} : two \emph{list} containing the possibles coordinates in each dimension.
}
}
\examples{
#### 1- with simulated data ####
## simulate
set.seed(10)
n <- 4
Y <- rnorm(n^2)

## conversion
res1 <- df2array(contrast = Y, coords = expand.grid(1:n + 0.5, 1:n + 0.5))
res2 <- df2array(contrast = Y, coords = expand.grid(1:n, 1:n), format = "matrix")
res3 <- df2array(contrast = Y, coords = expand.grid(2 * (1:n), 2 * (1:n)))
res4 <- df2array(contrast=cbind(Y ,Y, Y), coords = expand.grid(2 * (1:n), 2 * (1:n)),
                    range.coords = c(10,10))
					
## display
par(mfrow = c(2,2), mar = rep(2,4), mgp = c(1.5,0.5,0))
fields::image.plot(unique(res1$coords[,1]), unique(res1$coords[,2]), res1$contrast[[1]],
                   xlab = "", ylab = "")
fields::image.plot(unique(res2$coords[,1]), unique(res2$coords[,2]), res2$contrast,
                   xlab = "", ylab = "")
fields::image.plot(res3$contrast[[1]])
fields::image.plot(res4$contrast[[2]])

#### 2- with MRIaggr data ####
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")
carto <- selectContrast(MRIaggr.Pat1_red, param = "DWI_t0", format = "vector")
coords <- selectCoords(MRIaggr.Pat1_red)
coords[,1] <- coords[,1] + 30
coords[,2] <- coords[,2] + 15

## converion 1
array.DWI_t0 <- df2array(carto, coords = coords, default_value = 1000)$contrast[[1]]

# display
fields::image.plot(min(coords[,1]):max(coords[,1]), min(coords[,2]):max(coords[,2]),
                   array.DWI_t0[,,1], xlab = "i", ylab = "j")

## conversion 2
array.DWI_t0 <- df2array(contrast=carto, coords = coords, default_value = 1000, 
                         range.coords = c(128,128,3))$contrast[[1]]

# display
fields::image.plot(1:128, 1:128, array.DWI_t0[,,1], xlab = "i", ylab = "k")

}
\concept{convert.}
\keyword{functions}
