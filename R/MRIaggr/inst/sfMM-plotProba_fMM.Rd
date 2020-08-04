\name{plotProba_fMM}
\title{Graphical display of the probability distribution by group}
\alias{plotProba_fMM}
\usage{
plotProba_fMM(res.EM, type="posterior_proba",
    main="group", col=c("red","blue","green","yellow","purple","brown"),
    border="black", cex.main=0.75, window=FALSE,
    filename="EM.plotProb", width=480, height=480, path=NULL, unit="px", res=NA)
}
\description{
Graphical display of the probability distribution by group using the result of \code{\link{fMMalgo}} algorithm.
}
\arguments{
  \item{res.EM}{the result of the \code{\link{fMMalgo}} or the \code{\link{launcher_fMM}} function.}
  \item{type}{the type of proba to be display : can be \code{"prior_proba"} or \code{"posterior_proba"}}
  \item{main}{a main title for the plot (same as in \code{\link{plot.default}}). Default is \code{ group}.}
  \item{col}{the colors for lines (same as in \code{\link{plot.default}}). Default is \code{c("red","blue","green","yellow","purple","brown")}.}
  \item{border}{the color of the border around the bars (same as in \code{\link{hist}}). Default is \code{"black"}.}
  \item{cex.main}{the magnification to be used for main titles relative to the current setting of \code{cex}. (same as in \code{\link{plot.default}}).}
  \item{window}{Display the graphics in the current graphics device (\code{F}), in a new device (\code{T}), or save it in a png/eps/svg file (\code{"png"} / \code{"esp"} / \code{"png"}). \code{NULL} indicates no display.}
  \item{filename}{Name of the file (without the extension) where the plot should be saved. \code{character}. Ignored if window is not png/eps/svg.}
  \item{width}{the width of the device. \code{numeric}. Ignored if window is not png/eps/svg.}
  \item{height}{the width of the device. \code{numeric}. Ignored if window is not png/eps/svg.}
  \item{path}{the directory where the plot file will be created. \emph{character}. Default is \code{NULL} corresponding to the current working directory.}
  \item{unit}{the units in which \code{height} and \code{width} are given. Can be \code{"px"}, \code{"in"}, \code{"cm"} or \code{"mm"}. Default is \code{"px"}.}
  \item{res}{the nominal resolution in ppi which will be recorded in the bitmap file, if a positive integer. \emph{numeric}. Default is \code{NA}.}
}
                                                                                                                                      
\examples{
data(Simulation_data,package="fMMseg")

res_EMbase <- launcher_fMM(G=3,data=Simulation_data,
                          var_reg="Y",family=gaussian(link="identity"),
                          epsilon=5*10^{-3})

plotProba_fMM(res_EMbase)
}



                        
                        
                      