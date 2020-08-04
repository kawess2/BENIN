\name{initParameter}
\title{parameter initialization}
\alias{initParameter}
\alias{initParameter,MRIaggr-method}
\description{
Check and initialize the \code{param} argument in MRIaggr methods. For internal use.
}
\usage{
\S4method{initParameter}{MRIaggr}(object, param, test = TRUE, init = FALSE,
         accept.coords = TRUE, accept.mask = TRUE, accept.index = TRUE,
         arg_name = "param", long_name = "parameters", method)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}.}
  \item{param}{the contrast parameters to check or initialize. \emph{character vector} or \code{NULL}.}
  \item{test}{should the parameters be checked ? \emph{logical}.}
  \item{init}{should the parameters be initialized if \code{param} equals \code{NULL} ? \emph{logical}.}
  \item{accept.coords}{should coordinates be accepted as parameters ? \emph{logical}.}
  \item{accept.mask}{should mask be accepted as a parameter ? \emph{logical}.}
  \item{accept.index}{should index be accepted as a parameter ? \emph{logical}.}
  \item{arg_name}{a short name for the error message. \emph{character}.}
  \item{long_name}{the complete name for the error message. \emph{character}.}
  \item{method}{the name of the function that called the initializer. \emph{character}.}
}
\details{
ARGUMENTS : \cr
Setting \code{param} to \code{NULL} leads to load all available parameters including coodinates, mask and index if \code{accept.coords}, \code{accept.mask} or \code{accept.index} are respectively set to \code{TRUE}.
}
\concept{init.}
\keyword{methods,internal}
