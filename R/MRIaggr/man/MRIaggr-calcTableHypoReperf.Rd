\name{calcTableHypoReperf}
\title{Compute reperfusion and hypoperfusion tables}
\alias{calcTableHypoReperf}
\alias{calcTableHypoReperf,MRIaggr-method}
\description{
Compute of the reperfusion and hypoperfusion values.
}
\usage{
\S4method{calcTableHypoReperf}{MRIaggr}(object, param, timepoint, threshold = 1:10, sep = "_",
     norm_mu = FALSE, norm_sigma = FALSE, mask = NULL, numeric2logical = FALSE,
     param.update = "reperf",
	 verbose = optionsMRIaggr("verbose"), update.object = FALSE, overwrite = FALSE)
}
\arguments{
  \item{object}{an object of class \code{\linkS4class{MRIaggr}}. REQUIRED.}
  \item{param}{the perfusion parameter(s). \emph{character vector}. REQUIRED.}
  \item{timepoint}{one or two time points. \emph{character vector}. REQUIRED.}
  \item{threshold}{the value of the hypoperfusion thresholds. \emph{numeric vector}.}
  \item{sep}{the separator between the parameter names and the time points. \emph{character}.}
 
  \item{mask}{the binary contrast parameter indentifing the lesion at timepoint[1]. \emph{character} or \code{NULL} if no mask is available.}
  \item{numeric2logical}{should \code{mask} be converted to logical ? \emph{logical}.}
  \item{norm_mu}{the type of centering to apply on the parameter values. \emph{character}.}
  \item{norm_sigma}{the type of scaling to apply on the parameter values. \emph{character}.}

  \item{verbose}{should the execution of the function be traced ? \emph{logical}.}
  \item{param.update}{which type of parameter should be stored in the object ? Any of \code{"shift"} \code{"reperf"} \code{"reperf_pc"} \code{"deperf"} \code{"deperf_pc"}.}
  \item{update.object}{should the resulting values be stored in \code{object} ? \emph{logical}.}
  \item{overwrite}{if reperfusion or hypoperfusion values are already stored in \code{object@table_reperfusion} or \code{object@table_hypofusion}, can they be overwritten ? \emph{logical}.}  
}
\details{
ARGUMENTS: \cr
Information about the \code{norm_mu} and \code{norm_sigma} arguments can be found in the details section of \code{\link{selectContrast}}.

FUNCTION: \cr
If \code{mask} is set to \code{NULL}, no mismatch volume will not be computed.
}
\seealso{
\code{\link{calcThresholdMRIaggr}} to process contrast parameters. \cr
\code{\link{selectTable}} to select the reperfusion/hypoperfusion tables. \cr
\code{\link{calcW}} to compute the neighboring matrix.
}
\value{
  An \emph{list} containing :
 \itemize{
 \item\code{[[voxel]]} : a \emph{data.frame} containing the coordinates and the reperfusion values.
 \item\code{[[volume_hypo]]} : the number of hypoperfused observations for the various thresholds.
 \item\code{[[volume_reperf]]} : the number of reperfused observations for the various thresholds.
 }
}
\examples{
## load a MRIaggr object
data("MRIaggr.Pat1_red", package = "MRIaggr")

#### 1- directly ####
res <- calcTableHypoReperf(MRIaggr.Pat1_red, param = c("TTP","MTT"), timepoint=c("t0","t1"),
     mask="MASK_DWI_t0", numeric2logical = TRUE,
     update.object = TRUE, overwrite = TRUE)

carto_TTP_t0 <- selectContrast(MRIaggr.Pat1_red, param = "TTP_t0")
carto_TTP_t1 <- selectContrast(MRIaggr.Pat1_red, param = "TTP_t1")

## hypoperfusion
sum( (carto_TTP_t0 >= 4) )
selectTable(MRIaggr.Pat1_red, "hypoperfusion")["4","Vhypo.TTP_t0"]

## mismatch
testN <- (selectContrast(MRIaggr.Pat1_red, param = "MASK_DWI_t0") == 0)

sum( (carto_TTP_t0 >= 4) * testN )
selectTable(MRIaggr.Pat1_red, "hypoperfusion")["4","Vmismatch.TTP"]

sum( (carto_TTP_t0 >= 4) * testN ) / sum( testN == FALSE )
selectTable(MRIaggr.Pat1_red, "hypoperfusion")["4","PCmismatch.TTP"]


## reperfusion
sum((carto_TTP_t0 >= 4) * (carto_TTP_t1 < 4))
selectTable(MRIaggr.Pat1_red, "reperfusion")["4","Vreperf.TTP"]

sum((carto_TTP_t0 >= 4) * (carto_TTP_t1 < 4)) / sum( (carto_TTP_t0>=4) )
selectTable(MRIaggr.Pat1_red, "reperfusion")["4","PCreperf.TTP"]

## W reperfusion
carto_TTPth_t0 <- carto_TTP_t0
carto_TTPth_t0[carto_TTPth_t0 > 10] <- 10
carto_TTPth_t0[carto_TTPth_t0 < 0] <- 0

carto_TTPth_t1 <- carto_TTP_t1
carto_TTPth_t1[carto_TTP_t1 > 10] <- 10
carto_TTPth_t1[carto_TTP_t1 < 0] <- 0

weight <- (carto_TTPth_t0 - carto_TTPth_t1) / carto_TTPth_t0
weight[ ((carto_TTPth_t0 == 0) + (carto_TTP_t0 < 4) + (carto_TTP_t1 >= 4)) > 0 ] <- 0

sum((carto_TTP_t0 >= 4) * (carto_TTP_t1 < 4) * weight)
selectTable(MRIaggr.Pat1_red, "reperfusion")["4","VreperfW.TTP"]

sum((carto_TTP_t0 >= 4) * (carto_TTP_t1 < 4) * weight) / sum( (carto_TTP_t0 >= 4) )
selectTable(MRIaggr.Pat1_red, "reperfusion")["4","PCreperfW.TTP"]


## deperfusion
sum((carto_TTP_t0 < 4) * (carto_TTP_t1 >= 4))
selectTable(MRIaggr.Pat1_red, "reperfusion")["4","Vdeperf.TTP"]

sum((carto_TTP_t0 < 4) * (carto_TTP_t1 >= 4)) / sum( (carto_TTP_t0 >= 4) )
selectTable(MRIaggr.Pat1_red, "reperfusion")["4","PCdeperf.TTP"]

## shift
sum((carto_TTPth_t0 - carto_TTPth_t1 >= 4))
selectTable(MRIaggr.Pat1_red, "reperfusion")["4","Vshift_reperf.TTP"]

sum((carto_TTPth_t0-carto_TTPth_t1 >= 4)) / sum( (carto_TTP_t0>=4) )
selectTable(MRIaggr.Pat1_red, "reperfusion")["4","PCshift_reperf.TTP"]

#### 2- via calcThresholdMRIaggr ####
## perform segmentation  (call mritc)
\dontrun{
calcThresholdMRIaggr(MRIaggr.Pat1_red, param = c("TTP_t0","MTT_t0","TTP_t1","MTT_t1"),
         threshold=1:10, name_newparam = c("TTP.GR_t0","MTT.GR_t0","TTP.GR_t1","MTT.GR_t1"),
         rm.CSF=TRUE, hemisphere = "lesion",
         GRalgo=TRUE, seed = c("MASK_T2_FLAIR_t2","MASK_DWI_t0"), W.range  = sqrt(2),
         update.object = TRUE, overwrite = TRUE)

res <- calcTableHypoReperf(MRIaggr.Pat1_red, param = c("TTP.GR","MTT.GR"),
         timepoint = c("t0","t1"), mask = "MASK_DWI_t0", numeric2logical = TRUE,
         update.object = TRUE, overwrite = TRUE)
		 
## display
selectTable(MRIaggr.Pat1_red, "hypoperfusion")["4","Vhypo.TTP.GR_t0"]

par(mfrow = c(2,4), mar = rep(1.5,4), mgp = c(2,0.5,0))
multiplot(MRIaggr.Pat1_red, param = "TTP_t0", num = 1:3,
             palette = rainbow(10), window = NULL, main = "raw - slice ", 
			 breaks=(0:10)-10^{-10})
multiplot(MRIaggr.Pat1_red, param = "TTP.GR_t0", num = 1:3,
             palette = rainbow(10), window = NULL, main = "GR - slice ", 
			 breaks=(0:10)-10^{-10})
}
}
\concept{calc.}
\keyword{methods}

