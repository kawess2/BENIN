\name{EDK}
\title{Gaussian kernel}
\alias{EDK}
\description{
Apply a gaussian kernel to observations. For internal use.
}
\usage{
EDK(x, bandwidth, power = 2)
}
\arguments{
  \item{x}{the data on which the kernel will be applied. \emph{numeric} or \emph{numeric vector}.}
  \item{bandwidth}{the bandwidth of the kernel. \emph{numeric}.}
  \item{power}{the power of the kernel. \emph{numeric}.}
}
\details{
FUNCTION: \cr
This function corresponds to the following kernel :
  \deqn{\frac{1}{\sqrt[power]{2*pi*bandwidth^2}} *\exp(-(x/bandwidth)^power)}
Setting \code{power} to \code{2} lead to use a gaussian kernel.
}
\value{
  A \emph{numeric}.
}
\keyword{functions}
