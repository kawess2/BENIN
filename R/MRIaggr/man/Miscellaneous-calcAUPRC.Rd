\name{calcAUPRC}
\title{Area under the PR curve}
\alias{calcAUPRC}
\description{
Compute the area under the precision recall curve by numerical integration.
}
\usage{
calcAUPRC(x, y, subdivisions = 10000, performance = NULL, ci = TRUE, alpha = 0.05, 
          method = "Kronrod", reltol = .Machine$double.eps^0.25)
}
\arguments{
  \item{x}{the biomarker values. \emph{numeric vector}. REQUIRED.}
  \item{y}{the class labels. \emph{numeric vector}, \emph{character vector} or \emph{logical vector}. REQUIRED.}
  \item{subdivisions}{the maximum number of subintervals used for the integration. \emph{positive integer}. Only used if \code{method="integrate"}.}
  \item{performance}{an object of class \code{performance} can be supplied instead of arguments \code{x} and \code{y}.}  
  \item{ci}{should the confidence interval be computed ? \emph{logical}.}  
  \item{alpha}{the type 1 error rate. \emph{numeric}.}  
  \item{method}{the integration method used to compute the area under the curve. Any of \code{"integrate"}, \code{"Kronrod"}, \code{"Richardson"} \code{"Clenshaw"}, \code{"Simpson"} or \code{"Romberg"}.}  
  \item{reltol}{the relative accuracy requested. Positive \emph{numeric}.}  
}
\details{
This function requires to have installed the \emph{ROCR} package to work. \cr
The numeric integration of the precision over the recall values can be performed either using the \emph{integrate} function of the stats package (if \code{method="integrate"}) or using the \emph{integral} function of the pracma package. In the latter case, the method argument is used to define the integration procedure (see the documentation of \emph{integral} for more details). \cr
The confidence interval is computed using the first order delta method and the logistic transformation :
\deqn{ IC(AUPRC) = \left[ \frac{e^{\mu_\eta - 1.96 \tau}}{1+ e^{\mu_\eta - 1.96 \tau}} \; ; \; \frac{e^{\mu_\eta + 1.96 \tau}}{1+ e^{\mu_\eta + 1.96 \tau}} \right] }
\deqn{\mu_\eta = logit(\widehat{AUPRC})}
\deqn{\tau = \frac{1}{\sqrt{n*\widehat{AUPRC}*(1-\widehat{AUPRC})}}}
See section 3.2 of (Boyd, 2013) for more details.

ARGUMENTS: \cr
\code{y} must have exactly two levels.

If \code{performance} is set to \code{NULL}, the code{x} and \code{y} will be used to form the \code{performance} object.
}
\references{
Kendrick Boyd1, Kevin H. Eng, and C. David Page. \emph{Area Under the Precision-Recall Curve: Point Estimates and Confidence Intervals}. Machine Learning and Knowledge Discovery in Databases, 2013:451-466.
}
\value{
  If \code{ci=FALSE} a \emph{numeric} between 0 and 1. \cr
  If \code{ci=TRUE} a \emph{numeric vector} of length 3 containing the punctual estimation, the lower and the upper bound of the confidence interval.
}
\examples{
data(MRIaggr.Pat1_red, package = "MRIaggr")

## select parameter and binary outcome
cartoT2 <- selectContrast(MRIaggr.Pat1_red, param = "T2_FLAIR_t2", format = "vector")
cartoMASK <- selectContrast(MRIaggr.Pat1_red, param = "MASK_T2_FLAIR_t2", format = "vector")

## compute AUPRC
T2.AUPRC <- calcAUPRC(x = cartoT2, y = cartoMASK)

## compute AUC
\dontrun{
if(require(pROC)){
T2.AUC <- auc(roc(cartoMASK ~ cartoT2))
} 


## display
multiplot(MRIaggr.Pat1_red,param = "T2_FLAIR_t2", num = 1,
          index1 = list(coords = "MASK_T2_FLAIR_t2", outline = TRUE)
)
}

#### 2- with simulated data ####
n0 <- 1000
n1 <- c(10,100,1000)
for(iter_n in 1:length(n1)){
  X <- c(rnorm(n0,0),rnorm(n1[iter_n],2))
  Y <- c(rep(0,n0),rep(1,n1[iter_n]))
  print(calcAUPRC(X,Y))
}

## alternative way using a performance object
perfXY <- ROCR::performance(ROCR::prediction(X,Y), x.measure = "rec", measure = "prec")
calcAUPRC(performance = perfXY, subdivisions = 10000)
}
\concept{const.}
\keyword{functions}

