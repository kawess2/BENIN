\name{plotCv_fMM}
\title{Graphical display of the convergence criteria}
\alias{plotCv_fMM}
\description{
Graphical display of the convergence criteria using the result of \code{\link{fMMalgo}} algorithm.
}
\usage{
plotCv_fMM(res.EM, window=FALSE,
    filename="EM.traceLv", width=480, height=480, path=NULL, unit="px", res=NA)
}
\arguments{
  \item{res.EM}{The result of the \code{\link{fMMalgo}} or the \code{\link{launcher_fMM}} function.}
  \item{window}{Display the graphics in the current graphics device (\code{F}), in a new device (\code{T}), or save it in a png/eps/svg file (\code{"png"} / \code{"esp"} / \code{"png"}). \code{NULL} indicates no display.}
  \item{filename}{Name of the file (without the extension) where the plot should be saved. \code{character}. Ignored if window is not png/eps/svg.}
  \item{width}{the width of the device. \code{numeric}. Ignored if window is not png/eps/svg.}
  \item{height}{the width of the device. \code{numeric}. Ignored if window is not png/eps/svg.}
  \item{path}{the directory where the plot file will be created. \emph{character}. Default is \code{NULL} corresponding to the current working directory.}
  \item{unit}{the units in which \code{height} and \code{width} are given. Can be \code{"px"}, \code{"in"}, \code{"cm"} or \code{"mm"}. Default is \code{"px"}.}
  \item{res}{the nominal resolution in ppi which will be recorded in the bitmap file, if a positive integer. \emph{numeric}. Default is \code{NA}.}
}
\value{
  A \code{list} containing :
  \item{diff_lv}{a \code{vector} indicating the evolution of the likelihood between two consecutive steps.}
  \item{test.croissance}{a \code{logical} testing the monotonicity of the likelihood across steps.}
  \item{diff_lv_completee}{a \code{vector} indicating the evolution of the completed likelihood between two consecutive steps.}
  \item{test.croissance_completee}{a \code{logical} testing the monotonicity of the completed likelihood across steps.}
}
\examples{
data(Simulation_data,package="fMMseg")

res_EMbase <- launcher_fMM(G=3,data=Simulation_data,
                          var_reg="Y",family=gaussian(link="identity"),
                          epsilon=5*10^{-3})

plotCv_fMM(res_EMbase)
}



                        
                        
                      