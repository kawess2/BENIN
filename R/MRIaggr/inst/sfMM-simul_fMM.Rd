\name{simul_fMM}
\title{Simulation under the mixture model}
\alias{simul_fMM}
\usage{
simul_fMM(n, G, Ydist="gaussian", Y_mu, Y_sd, nclass = "auto",
    beta=FALSE, Xdist=NULL, X_mu=NULL, X_sd=NULL,
    n.Id=NULL, Id_sd=NULL, noise_sd=NULL,
    display_Y=TRUE, display_X=TRUE, window=FALSE,
	filename="calc-echantillon", width=1000, height=700, path=NULL, unit="px", res=NA)
}						
\description{
Simulation of a sample following a distribution issue from a non-spatial mixture model.
}
\arguments{
  \item{n}{The number of observation to simulate. \code{interger}.}
  \item{G}{The number of groups to be considered for the mixture. \code{integer}.}
  \item{Ydist}{The distribution of the observations. \code{character vector[G]}. Must contains only : \code{"unif"}, \code{"unif_segment"}, \code{"gaussian"}, \code{"Gamma"} or \code{"Beta"}.}
  \item{Y_mu}{The expected value of each group. \code{numeric vector[G]}.}
  \item{Y_sd}{The standard deviation of each group. \code{numeric vector[G]}.}
  \item{nclass}{The number of observation by class for the histogram. \code{integer}.}
  \item{beta}{Should a concomitant model be considered ? \code{logical}.}
  \item{Xdist}{The distribution of the concomitant variable. \code{character vector[G]}. Must contains only : \code{"unif"}, \code{"gaussian"}, \code{"Gamma"}.}
  \item{X_mu}{The expected value of the concomitant varible for each group. \code{numeric vector[G]}.}
  \item{n.Id}{The number of patient to be considered. \code{numeric}.}
  \item{Id_sd}{The standard deviation of the random effect. \code{numeric}.}
  \item{X_sd}{The standard deviation of the concomitant varible for each group. \code{numeric vector[G]}.}
  \item{noise_sd}{The standard deviation of the noise. \code{numeric}. \code{NULL} indicates no noise.}
  \item{display_Y}{Should the graphics representing the Y distribution be displayed. \code{logical}.}
  \item{display_X}{Should the graphics representing the X distribution be displayed. \code{logical}.}
  \item{window}{Display the graphics in the current graphics device (\code{F}), in a new device (\code{T}), save it in a png/eps/svg file (\code{"png"} / \code{"esp"} / \code{"png"}). \code{NULL} indicates no display.}
  
  \item{filename}{Name of the file (without the extension) where the plot should be saved. \code{character}. Ignored if window is not png/eps/svg.}
  \item{width}{the width of the device. \code{numeric}. Ignored if window is not png/eps/svg.}
  \item{height}{the width of the device. \code{numeric}. Ignored if window is not png/eps/svg.}
  \item{path}{the directory where the plot file will be created. \emph{character}. Default is \code{NULL} corresponding to the current working directory.}
  \item{unit}{the units in which \code{height} and \code{width} are given. Can be \code{"px"}, \code{"in"}, \code{"cm"} or \code{"mm"}. Default is \code{"px"}.}
  \item{res}{the nominal resolution in ppi which will be recorded in the bitmap file, if a positive integer. \emph{numeric}. Default is \code{NA}.}
}
\value{
   A \code{list} containing :
  \item{Y}{a \code{vector} containing the simulated response variables.}
  \item{X}{a \code{vector} containing the simulated concomitant variables.}
  \item{moy}{a \code{numeric} containing the empirical mean of the response variable.}
  \item{var}{a \code{numeric} containing the empirical standard deviation of the response variable.}
  \item{Y_g}{a \code{list} containing the simulated response variables for each group.}
  \item{X_g}{a \code{list} containing the simulated concomitant variables for each group.}
  \item{moy_g}{a \code{vector} containing the empirical mean of the response variable for each group.}
  \item{var_g}{a \code{vector} containing the empirical standard deviation of the response variable for each group.}
  \item{Ydist}{a \code{list} containing the family used to simulate the response variables.}
  \item{Y_mu}{a \code{vector} containing the simulated mean of the response variable for each group (without possible addition of noise).}
  \item{Y_sd}{a \code{vector} containing the simulated standard deviation of the response variable for each group (without possible addition of noise).}
  \item{beta}{a \code{vector} containing the theorical parameters of the concomitant model for each group.}
  \item{Xdist}{a \code{vector} containing the family used to simulate the concomitant variables.}
  \item{X_mu}{a \code{vector} containing the empirical expected value of the concomitant variable for each group.}
  \item{X_sd}{a \code{vector} containing the empirical standard deviation of the concomitant variable for each group.}
  \item{group}{a \code{matrix} indicating the group membership of the observations.}
  \item{prior}{a \code{matrix} indicating the prior probabilities associated to the concomitant model.}  
}
\examples{
# settings 
G <- 3 # nb of groups
n_px <- 30 # width of the area
n <- (n_px*G)^2 # total number of pixel
mu <- c(-3,0,4) # mean values of the groups
sigma <- c(sqrt(3),1,sqrt(3)) # sd values of the groups
Ydist <- c("gaussian","gaussian","gaussian") # distribution of the groups

# simulation 
set.seed(10)
simulation <- simul_fMM(n=rep(n/3,G),
                          G=G,                          
                          Ydist=Ydist,
                          Y_mu=mu,
                          Y_sd=sigma,
                          display_Y=TRUE)
                          
if(try(require(fields))){
image.plot(x=1:sqrt(n),y=1:sqrt(n),xlab="",ylab="",axes=FALSE,
           z=matrix(simulation$Y,sqrt(n),sqrt(n)),
           main="Groups")
}
}



                        
                        
                      