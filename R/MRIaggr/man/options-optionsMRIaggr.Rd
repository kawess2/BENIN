\name{optionsMRIaggr}
\title{Set or Query Default Values for MRIaggr}
\alias{optionsMRIaggr}
\description{
\code{optionsMRIaggr} works in the same way as \code{options} : it can be used to set or query the default values used by the functions of the MRIaggr package.
}
\usage{
optionsMRIaggr(..., reinit.options = FALSE)
}
\arguments{
  \item{...}{arguments in tag = value form, or a list of tagged values. \emph{character} or \emph{list}.}
  \item{reinit.options}{should the default values be reinitialized. \emph{logical}} 
}
\section{optionsMRIaggr elements}{

The individual optionsMRIaggr elements are:
\tabular{ll}{
  asp               \tab the aspect ratio y/x. \emph{numeric}. \cr
  axes              \tab should the axes be plotted ? \emph{logical}. \cr
  bg                \tab the color used for the background. \emph{character}. \cr
  breaks            \tab the break points or the number of breakpoints to use to generate the color intervals . \emph{numeric vector} or \emph{postive integer}. Active in \code{\link{multiplot}} only if argument \code{col} is \code{NULL}. \cr
  cex               \tab the expansion factor used to plot the observations. \emph{positive numeric}. \cr
  cex.index         \tab the expansion factor used to plot each additionnal points to display. \emph{integer vector of size 3}. \cr
  cex.legend        \tab the expansion factor of the legend. \emph{positive numeric}. \cr
  cex.main          \tab the expansion factor for the main title. \emph{numeric}. \cr
  checkArguments    \tab should the validity of thes argments be checked ? \emph{logical}. \cr
  col.index         \tab the color used to plot each additionnal points to display. \emph{integer vector of size 3}. \cr
  col.midplane      \tab the color used to display the midsagittal plan. \emph{character}. \cr 
  col.NA            \tab the color to use to plot the NAs. \emph{character}. \cr
  digit.legend      \tab the number of decimal places to use when displaying graphics. \emph{integer}. \cr
  digit.result      \tab the number of decimal places to use when displaying results. \emph{integer}. \cr
  digit.epsilon     \tab the number of decimal places to use when displaying convergence criteria. \emph{integer}. \cr
  digit.percentage  \tab the number of decimal places to use when displaying percentage. \emph{integer}. \cr
  filter.index      \tab the filter used to define the edge of a spatial group. Any of \code{"2D_N4"}, \code{"2D_N8"}, \code{"3D_N4"}, \code{"3D_N6"}, \code{"3D_N8"}, \code{"3D_N10"}, \code{"3D_N18"}, \code{"3D_N26"}. \cr
  height            \tab the height of the device used to export the plot. \emph{postive numeric}. \cr
  hemisphere        \tab the hemisphere to display. \emph{character}. \cr
  legend            \tab how should the legend be displayed ? \emph{logical}, \code{NULL} or \code{"only"}. \cr
  main              \tab an overall title for the plot. \emph{character}. \cr
  main.legend       \tab a main title for the legend. \emph{character}. \cr
  mar               \tab the number of margin lines to be specified on the four sides of the plot. \emph{positive numeric vector of size 4}. \cr
  mar.legend        \tab the number of margin lines to be specified on the four sides of the legend. \emph{numeric vector of size 4}. \cr
  mfrow             \tab the division of the device in plot region. \emph{numeric vector of size 2} or \code{NULL} leading automatic adjustment. \cr
  mgp               \tab the margin line for the axis title, axis labels and axis line. \emph{positive numeric vector of size 3}. \cr
  norm_mu           \tab the type of centering to apply on the parameter values. \emph{character}. \cr
  norm_sigma        \tab the type of scaling to apply on the parameter values. \emph{character}. \cr
  num.main          \tab should the slice number be written over each plot. \emph{logical}. \cr
  numeric2logical   \tab if a parameter is specified for the index arguments, should it be converted to logical ? \emph{logical}. \cr
  outline.index     \tab should all the additionnal points be displayed or only those on the boundary. \emph{logical}. \cr
  palette           \tab the colors or the palette to use when associating colors to observation intensities. \emph{character vector} or \emph{character}. Active in \code{\link{multiplot}} only if argument \code{col} is \code{NULL}. \cr
  path              \tab the directory where the plot file will be created. \emph{character}. If \code{path} is set to \code{NULL}, the image file is exported in the current working directory. \cr
  pch.index         \tab the label used to plot each additionnal points to display. \emph{integer vector of size 3}. \cr
  pch.NA            \tab the label to use to plot the NAs. \emph{postive integer}. \cr
  pty               \tab the type of plot region to be used. Can be \code{"s"} or \code{"m"}. \cr
  quantiles.legend  \tab should the quantiles of the data be displayed on the legend ? \emph{logical}. \cr
  res               \tab the nominal resolution in ppi which will be recorded in the bitmap file. \emph{positive integer}. \cr
  slice_var         \tab the type of view to use. \code{"i"} for sagittal view, \code{"j"} for coronal view and \code{"k"} for transverse view. \emph{character}. \cr
  type.breaks       \tab should the break points be equally space according the range of data values (\code{"range"}), centered (\code{"range_center"}) or correspond to the quantile values (\code{"quantile"}) ? \cr
  unit              \tab the units in which \code{height} and \code{width} are given. Can be any of \code{"px"}, \code{"in"}, \code{"cm"} or \code{"mm"}.. \cr
  verbose           \tab should the execution of the methods related to MRIaggr objects be traced ? \emph{logical}. \cr
  xlab              \tab a title for the x axis. \emph{character}. \cr
  ylab              \tab a title for the y axis. \emph{character}. \cr
  width             \tab the width of the device used to export the plot. \emph{postive numeric}. \cr
  window            \tab the type of device on which the plot will be displayed. \emph{logical}, \code{NULL} or \emph{character}. \cr
  }
}
\section{Details}{
ARGUMENTS: \cr
\code{palette} must have be a \emph{charater vector} of length \code{length(breaks)-1} containing colors or a single \emph{character} corresponding to a palette name : \code{"rainbow"}, \code{"grey.colors"}, \code{"heat.colors"}, \code{"terrain.colors"}, \code{"topo.colors"} or \code{"cm.colors"}.

Possible values for \code{legend} are:
\itemize{
\item \code{TRUE} : the legend is displayed in the current device.
\item \code{NULL} : the legend is displayed it in a new device.
\item \code{FALSE} : the legend is not displayed.
\item \code{"only"} : the legend is displayed but not the slices.
}

Possible values for \code{windows} are:
\itemize{
\item \code{NULL} : the plot is displayed on the current device with no reshape. 
\item \code{FALSE} : the plot is displayed on the current device with the appropriate reshape.
\item \code{TRUE} : a new graphical device is open.
\item \code{"eps"} : the plot is displayed in an image file using \code{\link{postscript}}.
\item \code{"svg"} : the plot is displayed in an image file using \code{\link{svg}}.
\item \code{"png"} : the plot is displayed in an image file using \code{\link{png}}.
\item \code{"pdf"} : the plot is displayed in an image file using \code{\link{pdf}}.
}

Arguments \code{filename}, \code{width}, \code{height}, \code{path}, \code{unit} and \code{res} correspond to those of \code{\link{postscript}}, \code{\link{svg}} and \code{\link{png}}.

Arguments \code{mfrow}, \code{mar}, \code{mar.legend}, \code{mgp} \code{pty}, \code{asp} and \code{bg} acorrespond to those of \code{\link{par}}.

See the details section of \code{\link{selectContrast}} for more details about the \code{hemisphere}, \code{norm_mu} and \code{norm_mu} arguments.

See the \code{\link{initNeighborhood}} function for details about \code{filter.index} possible values.

ARGUMENTS: \cr
The MRIaggr environment store an object called \code{.ls_optionsMRIaggr} that contains the list of default values. This object can be modified using the \code{optionsMRIaggr} function. The original default values can be restored using \code{optionsMRIaggr(reinit.options=TRUE)}.
}
\value{
If not argument is provided, a list containing all the graphical parameters is returned. \cr
If one argument is provided with no value, the value of the corresponding graphical parameter. \cr
If several arguments are provided, a list containing the corresponding graphical parameters is returned. \cr
If one argument is provided with a value, nothing is returned. \cr
}
\examples{
## load a MRIaggr object
data(MRIaggr.Pat1_red, package = "MRIaggr")

# display 3 slices 
multiplot(MRIaggr.Pat1_red, param = "DWI_t0",              
          num = 1:3, axes = FALSE, bg = "white", legend = FALSE)

##
# arguments :
# axes = FALSE, bg = "white", legend = FALSE 
# can become the default argument of multiplot using :
##

optionsMRIaggr(axes = FALSE, bg = "white", legend = FALSE)

multiplot(MRIaggr.Pat1_red, param = "DWI_t0")

# now axes can be restored using either
multiplot(MRIaggr.Pat1_red, param = "DWI_t0", axes = TRUE)
# or
optionsMRIaggr(axes = TRUE)
multiplot(MRIaggr.Pat1_red, param = "DWI_t0")

## to get default arguments of multiplot
optionsMRIaggr()
optionsMRIaggr()$bg
optionsMRIaggr("bg")

## to restore the original defaut values
optionsMRIaggr(reinit.options = TRUE)
}
\concept{plot.}
\keyword{functions}
