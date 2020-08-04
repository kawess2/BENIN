\name{tnorm}
\alias{dtnorm}
\alias{rtnorm}
\title{Truncated Normal distribution}
\description{Density and random generation for the truncated Normal distribution with mean equal to mean and standard deviation equal to sd before truncation, and truncated on the interval [lower, upper].
For internal use.}
\usage{
 dtnorm(x, mean = 0, sd = 1, lower = -Inf, upper = Inf, log = FALSE)
 rtnorm(n, mean = 0, sd = 1, lower = -Inf, upper = Inf)
}
\arguments{
\item{x}{vector of quantiles. \emph{numeric vector}.}
\item{n}{number of observations. \emph{integer}.}
\item{mean}{mean   vector of means. \emph{numeric}. Default is \code{0}.}
\item{sd}{vector of standard deviations. \emph{numeric}. Default is \code{1}.}
\item{lower}{truncation point. \emph{numeric}. Default is \code{-Inf}.}
\item{upper}{upper truncation point. \emph{numeric}. Default is \code{Inf}.}
\item{log}{if TRUE, probabilities p are given as log(p). \emph{logical}. Default is \code{FALSE}.}
}
\details{
These functions are a copy of the \code{dtnorm} and \code{rtnorm} of the \emph{msm} package. See the corresponding help page for additional details.
}
