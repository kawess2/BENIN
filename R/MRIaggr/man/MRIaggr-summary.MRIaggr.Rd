\name{summary-methods}
\title{Summary Method for Class "MRIaggr"}
\alias{show,MRIaggr-method}
\alias{summary,MRIaggr-method}
\description{
 Summarize some information of an object of class "MRIaggr". 
}
\usage{
\S4method{show}{MRIaggr}(object)

\S4method{summary}{MRIaggr}(object, param = FALSE, clinic = FALSE, 
         descStats = FALSE, history = FALSE, verbose = optionsMRIaggr("verbose"))
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{param}{should detailed information be printed for the contrast parameters ? \emph{logical}.}
  \item{clinic}{should detailed information be printed for the clinical attribute ? \emph{logical}.}
  \item{descStats}{should detailed information be printed for the \code{ls_descStats} attribute ? \emph{logical}.}
  \item{history}{should the\code{calc} and \code{const} methods that have been applied to the object be listed ? \emph{logical}.}
  \item{verbose}{should a summary of all attributes should be displayed or only those which are non-empty ? \emph{logical}.}
}
\value{
  None.
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")
show(MRIaggr.Pat1_red)

# display all elements
summary(MRIaggr.Pat1_red)

# display non-empty elements
summary(MRIaggr.Pat1_red, verbose = FALSE)

# display a summary of the contrast parameters
summary(MRIaggr.Pat1_red, param = TRUE)

# display a summary of the clinical data
summary(MRIaggr.Pat1_red, clinic = TRUE)

# display a summary of the elements in descStats
summary(MRIaggr.Pat1_red, descStats = TRUE)

# display processing methods previously applied on the object
summary(MRIaggr.Pat1_red, history = TRUE)

}
\concept{plot.}
\keyword{methods,summary}

