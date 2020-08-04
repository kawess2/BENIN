\name{plotY_fMM}
\alias{plotY_fMM}
\title{Graphical display of the mixture distribution}
\description{Graphical display of the mixture distribution using the result of \code{\link{fMMalgo}} algorithm.}
\usage{
plotY_fMM(res.EM, prob="posterior_proba", type="density", numY=1, identifier=NULL, nclass=NULL, 
    x.legend="topright", col=NULL, main=NULL, cex.main=1, lwd=1, lwd.mean=2, lty.wmean=2, 
    window=FALSE, filename="plotY_fMM", width=480, height=480, path=NULL, unit="px", res=NA)
}
\arguments{
  \item{res.EM}{The result of the \code{\link{fMMalgo}} or the \code{\link{launcher_fMM}} function.}
  \item{prob}{the type of group membership : can be \code{"group"}, \code{"prior_proba"} or \code{"posterior_proba"}.}
  \item{type}{the type of plot : can be \code{"hist"} (only available if \code{prob="group"}) or \code{"density"}.}
  \item{numY}{the number of the response to be selected. \code{integer}. Default is \code{1}.}
  \item{identifier}{the identifier of observation to be displayed. \code{NULL} indicates all observations.}
  \item{nclass}{a single number giving the number of cells for the histogram (same as in \code{\link{hist}}).}
  \item{x.legend}{the position of the legend by keywords (same as the \code{x} argument of \code{\link{legend}}), possible values are : "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center".}
  \item{col}{the colors for lines (same as in \code{\link{plot.default}}). Default is \code{c("red","blue","green","yellow","purple","brown")}.}
  \item{main}{a main title for the plot (same as in \code{\link{plot.default}}). Default is \code{ group}.}
  \item{cex.main}{the magnification to be used for main titles relative to the current setting of \code{cex}. (same as in \code{\link{plot.default}}).}
  \item{lwd}{a vector of line widths for the density lines (same as in \code{\link{plot.default}}).}
  \item{lwd.mean}{a vector of line widths for the mean lines (same as in \code{\link{plot.default}}).}
  \item{lty.wmean}{a vector of line types for the mean lines (same as in \code{\link{plot.default}}).}
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
                          epsilon=5*10^{-3},iter_max=100,
                          trace_iter=FALSE)

plotY_fMM(res_EMbase)
}
