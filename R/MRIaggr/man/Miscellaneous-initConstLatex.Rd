\name{initConstLatex}
\title{Initializateurs for the constLatex function}
\alias{initDirPat_constLatex}
\alias{initPlot_constLatex}
\alias{initSection_constLatex}
\usage{
initDirPat_constLatex(directory, param, identifier, verbose)

initSection_constLatex(directory.plot, names_dirs, param, subsection, label)

initPlot_constLatex(directory, names_dirs, directory.plot, identifier, 
         tabular, plotPerPage)

}
\description{
Initialize arguments of the constLatex function. For internal use.
}
\arguments{
  \item{directory}{the path to the root directory. This directory should contains subsdirectories themself containing the image files. \emph{character}.}
  \item{directory.plot}{the paths to the subdirectory(ies) containing the image files. \emph{character} or \emph{character vector}.}
  \item{names_dirs}{the names of the subdirectory(ies). \emph{character} or \emph{character vector}.}
  \item{identifier}{the identifiers of the patients for which the graphics should be displayed. \emph{character vector} or \code{NULL} leading to use all patients.}
  \item{param}{the names of directories containing the images. \emph{character vector} or \code{NULL} leading to use all directories.}
  \item{tabular}{a list of data.frame to display in the table format. \emph{list of data.frame} or \code{NULL} if there is no table to display.}
  \item{subsection}{the names of subsections for the latex document. \emph{character vector} or \code{NULL} leading to use \code{param} for naming the subsections.}
  \item{label}{the legend that will be displayed under each figure of a given parameter. \emph{character vector} or \code{NULL} leading to use \code{param} for the legend.}
  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
  \item{plotPerPage}{the number of image that should be displayed on the same page. \emph{strictly positive interger}.}
  
}
\concept{init.}
\keyword{functions,internal}
