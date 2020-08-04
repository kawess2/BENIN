\name{calcCriteriaGR}
\title{Assessment of clustering quality}
\alias{calcCriteriaGR}
\description{
Compute several quality indexes of a two group clustering. For internal use.
}
\usage{
calcCriteriaGR(contrast, groups, W = NULL, sigma = NULL, breaks = NULL, 
         rm.warning = TRUE, criterion.transition = FALSE, criterion.sdfront = FALSE,
         criterion.entropy = TRUE, criterion.Kalinsky = TRUE, criterion.Laboure = TRUE)
}
\arguments{
  \item{contrast}{the contrast value of each observations. \emph{numeric vector}. REQUIRED.}
  \item{groups}{the indicator of group membership. \emph{logical vector}. REQUIRED.}
  \item{W}{the neighbourhood matrix. \emph{dgCMatrix} or \code{NULL} leading to not compute the d1 criterion.}
  \item{sigma}{the \code{sigma_max} that have been used in the GR algorithm. \emph{positive numeric vector}.}
  \item{breaks}{the break points to use to categorize the contrast distribution. \emph{numeric vector}.}
  \item{rm.warning}{should warning be displayed. \emph{logical}.}
 
  \item{criterion.transition}{should the boundary criterion based on the transition levels be computed ? \emph{logical}. }  
  \item{criterion.sdfront}{should the boundary criterion based on the standard deviation be computed ? \emph{logical}.}  
  \item{criterion.entropy}{should the region criterion based on the entropy be computed ? \emph{logical}.}  
  \item{criterion.Kalinsky}{should the region criterion based on the Kalinsky index be computed ? \emph{logical}.}  
  \item{criterion.Laboure}{should the region criterion based on the Laboure index  be computed ? \emph{logical}.}
}
\references{
Chantal Revol-Muller, Francoise Peyrin, Yannick Carrillon and Christophe Odet. \emph{Automated 3D region growing algorithm based on an assessment function}. Pattern Recognition Letters, 23:137-150,2002.
}
\concept{calc.}
\keyword{functions,internal}
