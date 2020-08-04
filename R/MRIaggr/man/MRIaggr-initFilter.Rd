\name{initFilter}
\title{Initialization of a filter}
\alias{initFilter}
\usage{
initFilter(filter, method)
}
\description{
Return a filter corresponding to the specified name.
}
\arguments{
  \item{filter}{the filter to be initialized. \emph{character}. REQUIRED.}
  \item{method}{the name of the function that called the initializer. \emph{character}.}
}
\details{
ARGUMENTS : \cr
\code{filter} refers to classical filters that will be constructed by the function: \cr
The first two characters refer to the dimension of the filter : \code{"2D"} or \code{"3D"}. \cr
The third character must be \code{"_"}. \cr
The fourth character refers to the type of filter among : \code{"M"}, \code{"G"}, \code{"S"} or \code{"I"} :
 \itemize{
 \item \code{"M"} : median filter (by default a weight of 1 is attributed to each site)
 \item \code{"G"} : gaussian filter (sites are weight with a gaussian kernel)
 \item \code{"S"} : sobel filter (gradient in the i, j or k direction)
 \item \code{"I"} : influence filter (sites with distance to the central site inferior or equal to \code{sqrt(p)} have weight 1, otherwise 0)
 }
The last one or two characters indicates the size p of the filter except for the sobel filter where it indicates the direction of the gradient (\code{"x"}, \code{"y"} or \code{"z"})
}
\value{
A \emph{list} were : 
\itemize{
  \item\code{[[filter]]} : contains the filter. \emph{matrix of dimension p*p} or \emph{array of dimension p*p*p}.
  \item\code{[[filter_split]]} : contains the decomposition of the filter name. \emph{character vector}.
 }
}
\examples{
initFilter("2D_G3", method = "calcFilter")$filter
initFilter("2D_G5", method = "calcFilter")$filter
initFilter("3D_G3", method = "calcFilter")$filter

initFilter("2D_M9", method = "calcFilter")$filter
initFilter("3D_M3", method = "calcFilter")$filter

initFilter("2D_Sx", method = "calcFilter")$filter
initFilter("2D_Sy", method = "calcFilter")$filter
initFilter("3D_Sx", method = "calcFilter")$filter
initFilter("3D_Sz", method = "calcFilter")$filter

initFilter("2D_I3", method = "calcFilter")$filter
initFilter("2D_I5", method = "calcFilter")$filter
initFilter("2D_I7", method = "calcFilter")$filter
initFilter("2D_I9", method = "calcFilter")$filter
initFilter("3D_I5", method = "calcFilter")$filter
}
\concept{init.}
\keyword{functions}