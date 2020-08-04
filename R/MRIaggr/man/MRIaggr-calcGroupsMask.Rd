\name{calcGroupsMask}
\title{Compute spatial groups}
\alias{calcGroupsMask}
\alias{calcGroupsMask,MRIaggr-method}
\description{
Compute the spatial groups defined by a binary parameter.
}
\usage{
\S4method{calcGroupsMask}{MRIaggr}(object, mask, numeric2logical = FALSE,
     W = "ifany", W.range, W.spatial_res = c(1,1,1),
	 verbose = optionsMRIaggr("verbose"), update.object = FALSE, overwrite = TRUE)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{mask}{the binary contrast parameter that should be used to identifying the spatial groups. \emph{character}. REQUIRED.}
  \item{numeric2logical}{should \code{mask} be convert to logical ? \emph{logical}.}
  
  \item{W}{the neighbourhood matrix. \emph{dgCMatrix} or \code{"ifany"} leading to use the neighbourhood matrix stored in the \code{object} if any and else computate this matrix.}
  \item{W.range}{the neighbourhood range. \emph{postive numeric}. REQUIRED.}
  \item{W.spatial_res}{a dilatation factor for the coordinates. \emph{positive numeric vector of size 3}.}
  
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
  \item{update.object}{should the resulting spatial groups be stored in \code{object} ? \emph{logical}.}
  \item{overwrite}{if spatial groups are already stored in \code{object@ls_descStats}, can they be overwritten ? \emph{logical}.}  
}
\details{
This function requires to have installed the \emph{Matrix} and the \emph{spam} package to work.

FUNCTION: \cr
Call the \code{\link{calcGroupsW}} function.
}
\value{
  An \emph{list} containing for each parameter:
  \itemize{
  \item\code{[[group]]} : a \emph{vector} containing the group index for each observation.
  \item\code{[[group_size]]} : a \emph{vector} with the size of each spatial group.
  \item\code{[[group_number]]} : the number of spatial groups. \emph{integer vector}.
  \item\code{[[group_max]]} : the index of the largest group. \emph{integer vector}.
  }
}
\seealso{
\code{\link{selectDescStats}} to select the spatial groups. \cr
\code{\link{calcW}} to compute the neighbourhood matrix.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## compute spatial groups
calcGroupsMask(MRIaggr.Pat1_red, mask = c("MASK_DWI_t0","MASK_T2_FLAIR_t2"),
               W.range = 6, W.spatial_res = c(1.875,1.875,6),
               update.object = TRUE, overwrite = TRUE)

## extract spatial groups
selectDescStats(MRIaggr.Pat1_red, "GroupsLesion")
}
\concept{calc.}
\keyword{methods}

