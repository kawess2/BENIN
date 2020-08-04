\name{calcW}
\title{Compute the neighbourhood matrix}
\alias{calcW}
\alias{calcW,data.frame-method}
\alias{calcW,MRIaggr-method}
\description{
Computate a neighbourhood matrix using spatial coordinates.
}
\usage{
\S4method{calcW}{data.frame}(object, range, method = "euclidean", upper = NULL,
          format = "dgCMatrix", row.norm = FALSE, spatial_res = rep(1,ncol(object)),
		  calcBlockW = FALSE)

\S4method{calcW}{MRIaggr}(object, range, spatial_res = c(1,1,1), num = NULL,
     hemisphere = "both", subset = NULL, upper = TRUE, format = "dgCMatrix", 
	 row.norm = FALSE, calcBlockW = FALSE,
	 verbose = optionsMRIaggr("verbose"), update.object = FALSE, overwrite = FALSE)
}
\arguments{
  \item{object}{a \emph{data.frame} containing the coordinates of the observations or an \code{object} of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{range}{only distances smaller than \code{range} are stored in \code{W}. \code{postive numeric}. REQUIRED.}
  \item{method}{the distance measure to be used. \emph{character}. This must be one of \code{"euclidean"}, \code{"maximum"}, \code{"minkowski"} or \code{"greatcircle"}.}  
  \item{spatial_res}{a dilatation factor for the coordinates. \emph{positive numeric vector of size 3}.}
  \item{num}{the slices to use. \emph{numeric vector} or \code{NULL}.}
  \item{hemisphere}{the hemisphere to use. \emph{character}.}
  \item{subset}{the subset of observations to use. \emph{positive integer vector} or \code{NULL} leading to use all observations.}
  \item{upper}{should the entire matrix (\code{NULL}) or only the upper-triagonal (\code{TRUE}) or only the lower-triagonal (\code{FALSE}) values be calculated ?}
  \item{format}{the format of the output. Could be \code{"spam"} or \code{"dgCMatrix"}.}
  \item{row.norm}{should the resulting matrix be row-normalized ? \code{TRUE/FALSE}.}
  \item{calcBlockW}{should the partition into disjoint spatial blocks of sites be computed ? \emph{logical}.}
   
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
  \item{update.object}{should the resulting neighbourhood matrix be stored in \code{object} ? \emph{logical}.}
  \item{overwrite}{if a neighbourhood matrix is already stored in \code{object@ls_descStats}, can it be overwritten ? \emph{logical}.}
}
\details{
These functions require to have installed the \emph{Matrix} and the \emph{spam} package to work.

ARGUMENTS: \cr
Information about the \code{num} argument can be found in the details section of \code{\link{initNum}}.

The \code{range} argument cooresponds to the \code{delta} argument of the \code{nearest.dist} function from the \emph{spam} package.

Information about the \code{hemisphere} argument can be found in the details section of \code{\link{selectContrast}}.

The \code{row.norm} argument is ignored if \code{format} is set to \code{"spam"}.

FUNCTION: \cr
This function relies on the \code{nearest.dist} function of the \emph{spam} package.

Each of the \code{num}, \code{hemisphere} and \code{subset} argument define a subset of the total set of observations. 
It is the intersection of all these three subsets that is extracted.
}
\value{
  Invisible. A list containing :
  \itemize{
  \item\code{[[W]]} :  a \emph{spam} or \emph{dgCMatrix} object.
  \item\code{[[block]]} :  the output of the \code{calcBlockW} function if \code{calcBlockW} is set to \code{TRUE}, \code{NULL} otherwise.
}
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

#### 1- data.frame method ####
coords <- selectCoords(MRIaggr.Pat1_red, num = 1:3, hemisphere = "lesion")

## full W 
W <- calcW(object = coords, range = sqrt(2))$W
W[1:10,1:10]
table(spam::rowSums(W))

## full W normalized by row
W <- calcW(object = coords, range = sqrt(2), row.norm = TRUE)$W
W[1:10,1:10]
table(spam::rowSums(W))

## upper W 
W <- calcW(object = coords, range = sqrt(2), upper = TRUE)$W
W[1:10,1:10]

#### 2- MRIaggr method ####

## compute W (regular lattice)
W <- calcW(MRIaggr.Pat1_red, range = sqrt(2), upper = NULL, num = 1:3, hemisphere = "lesion")$W
table(spam::rowSums(W > 0))

## compute W (irregular lattice)
W <- calcW(MRIaggr.Pat1_red, range = sqrt(2 * 1.875^2), upper = NULL, num = 1:3, 
         hemisphere = "lesion", spatial_res=c(1.875, 1.875, 6))$W
table(spam::rowSums(W > 0))
}
\concept{calc.}
\keyword{methods}


                        
                        
                      
