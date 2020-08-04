\name{valid}
\title{Checking argument validity}
\alias{validCharacter}
\alias{validClass}
\alias{validDim_vector}
\alias{validDim_matrix}
\alias{validInteger}
\alias{validLogical}
\alias{validNames}
\alias{validNumeric}
\alias{validPath}

\description{
Check the validity of an argument. For internal use.
}
\usage{
validCharacter(value, name = as.character(substitute(value)), validLength, 
         validValues = "character", refuse.NULL = TRUE, method)

validClass(value, name = as.character(substitute(value)), validClass, 
         superClasses = TRUE, method)

validDim_vector(value1, value2, name1 = as.character(substitute(value1)), 
         name2 = as.character(substitute(value2)), type = "length", method)

validDim_matrix(value1, value2, name1 = as.character(substitute(value1)), 
         name2 = as.character(substitute(value2)), type = "both", method)

validInteger(value, name = as.character(substitute(value)), validLength, 
         validValues = NULL, min = NULL, max = NULL, 
         refuse.NA = TRUE, refuse.NULL = TRUE, refuse.duplicates = FALSE, method)
		 
validLogical(value, name = as.character(substitute(value)), validLength, 
         refuse.NULL = TRUE, refuse.NA = TRUE, method)		 

validNames(value, name = as.character(substitute(value)),
         validLength = NULL, validValues = NULL, method)

validNumeric(value, name = as.character(substitute(value)), validLength,
         validValues = NULL , min = NULL, max = NULL,
         refuse.NA = TRUE, refuse.NULL = TRUE, refuse.duplicates = FALSE, method)
		 
validPath(value, name = as.character(substitute(value)), method)
}
\arguments{
  \item{value, value1, value2}{the value of the argument(s).}
  \item{name, name1, name2}{the name(s) of the argument(s). \emph{character}.}
  \item{validLength}{the valid lengthes for the argument. \emph{integer vector}.}
  \item{validClass}{the valid classes for the argument. \emph{character vector}.}
  \item{validValues}{the valid values for the argument. \emph{character vector}.}
  \item{refuse.NA}{Is \code{NA} an invalid value ? \emph{logical}.}
  \item{refuse.NULL}{Is \code{NULL} an invalid value ? \emph{logical}.}
  \item{refuse.duplicates}{Are duplicated values invalid values. \emph{logical}.}
  \item{type}{the operator used to compare the size. \emph{character}. See the details section for more information.}
  \item{min}{the minimum value for the argument. \emph{numeric} or \emph{NULL}.}
  \item{max}{the maximum value for the argument. \emph{numeric} or \emph{NULL}.}
  \item{superClasses}{Should the super-classes of the object's class be also considered or only the object's class ? \emph{logical}.}
  \item{method}{the name of the function that called the valid function. \emph{character}.}
}
\details{
\code{validCharacter} checks whether \code{value} is a character vector.

\code{validClass} checks whether \code{value} belongs to the correct class of R objects.

\code{validDim_vector} checks whether \code{value1} is a vector matching the dimensions of \code{value2}. 
Argument \code{type} indicates which operator should be used to measure the size of \code{value2} among : \code{"length"}, \code{"nrow"} and \code{"ncol"}.

\code{validDim_matrix} checks whether \code{value1} is a matrix matching the dimensions of \code{value2}. 
Argument \code{type} indicates which operator should be used to measure the size of \code{value2} among : \code{"nrow"}, \code{"ncol"} and \code{"both"}, the latter indicating to use both \code{nrow} and \code{ncol}.

\code{validInteger} checks whether \code{value} is a integer vector.

\code{validLogical} checks whether \code{value} is a logical vector.

\code{validNames} checks whether the names \code{value} are valid.

\code{validNumeric} checks whether \code{value} is a numeric vector.

\code{validPath} checks whether \code{value} corresponds to a valid path.
}
\concept{valid.}
\keyword{function,internal}
