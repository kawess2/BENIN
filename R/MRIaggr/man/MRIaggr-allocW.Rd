\name{allocW}
\title{Allocate a neighbourhood matrix}
\alias{allocW<-}
\alias{allocW<-,MRIaggr-method}
\usage{
\S4method{allocW}{MRIaggr}(object, type, 
         overwrite = FALSE, verbose = optionsMRIaggr("verbose")) <- value
}
\description{
Allocate a neighbourhood matrix to a \code{\linkS4class{MRIaggr}} object.
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{value}{a list containing the neighbourhood matrix and associated elements. \emph{list}. See the details section. REQUIRED.}
  \item{type}{the type of element to be allocated. Can be \code{"Wmatrix"}, \code{"Wblocks"} or \code{"upper"}. REQUIRED.}
  \item{overwrite}{if an element with the same name is already stored in \code{object@ls_descStats}, can it be overwritten ? \emph{logical}.}  
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
}
\details{
ARGUMENTS: \cr
Argument \code{value} must be a list containing named elements :
\itemize{
\item Element \code{"Wmatrix"} must be matrix or a sparse matrix.
\item Element \code{"Wblocks"} must be consistent with the output of the \code{calcBlockW} method.
\item Element \code{"verbose"} must be logical or NULL method.
}

FUNCTION: \cr
This method enable to store a neighbourhood matrix (\code{type="Wmatrix"}), a block decomposition of the matrix (\code{type="Wblocks"})  and its format (\code{type="upper"}) : entire matrix (\code{NULL}), upper-triagonal (\code{TRUE}) or lower-triagonal (\code{FALSE}).
}
\seealso{
 \code{\link{calcW}} to compute the neighboring matrix. \cr
 \code{\link{calcBlockW}} to decompose the neighboring matrix in independant blocks. \cr
 \code{\link{selectW}} to extract the neighboring matrix or its related elements. 
}
\value{
  None.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## search for a neighbourhood matrix
# selectW(MRIaggr.Pat1_red) # no neighbourhood matrix

## compute the neighbourhood matrix
resW <- calcW(object = MRIaggr.Pat1_red, range = sqrt(2))

## store it
allocW(MRIaggr.Pat1_red, type = "Wmatrix") <- list(Wmatrix = resW$W)

## it is automatically performed when specifying update.object = TRUE
calcW(object = MRIaggr.Pat1_red, range = sqrt(2), update.object = TRUE, overwrite = TRUE)
}
\concept{alloc.}
\keyword{methods}
