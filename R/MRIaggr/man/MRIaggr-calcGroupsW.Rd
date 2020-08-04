\name{calcGroupsW}
\title{Compute spatial groups}
\alias{calcGroupsW}
\description{
Compute the spatial groups using a neighbourhood matrix.
}
\usage{
calcGroupsW(W, subset = NULL, max_groups = 10000, 
         verbose = optionsMRIaggr("verbose"))
}
\arguments{
  \item{W}{the neighbourhood matrix. \emph{dgCMatrix}. REQUIRED.}
  \item{subset}{the subset of observations to use. \emph{positive integer vector} or \code{NULL} leading to use all observations.}
  \item{max_groups}{the maximum number of groups. \emph{postive integer}.}
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
}
\details{
This function requires to have installed the \emph{Matrix} and the \emph{spam} package to work.
}
\value{
  An \emph{list} containing :
  \itemize{
  \item\code{[[group]]} : a \emph{vector} containing the group index for each observation. Observations out of the subset are set to \code{NA}.
  \item\code{[[subset]]} : a \emph{vector} containing the group index for each observation in subset.
  \item\code{[[group_size]]} : a \emph{vector} with the size of each spatial group.
  \item\code{[[group_number]]} : the number of spatial groups. \emph{integer vector}.
  \item\code{[[group_max]]} : the number of the largest group. \emph{integer vector}.
  }
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

## select data
MASK_DWI_t0 <- selectContrast(MRIaggr.Pat1_red, param = "MASK_DWI_t0")
coords <- selectCoords(MRIaggr.Pat1_red)

## select compute W
W <- calcW(object = as.data.frame(coords[MASK_DWI_t0 == 1,]), 
           range = sqrt(2), row.norm = TRUE)$W
 
## find spatial groups
res.Groups <- calcGroupsW(W)
res.Groups$group_size

## display
multiplot(coords[MASK_DWI_t0 == 1,], contrast = res.Groups$group,
          legend = FALSE, cex=0.5,
          palette=rainbow(10)[-1])
}
\concept{calc.}
\keyword{functions}

