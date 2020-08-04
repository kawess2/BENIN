\name{test_fMM}
\title{Checking of the \code{\link{fMMalgo}} arguments}
\alias{test_fMM}
\description{
Check the arguments of the \code{\link{fMMalgo}} function and return the appropriate error message if necessary.
For internal use.
}
\usage{
test_fMM(n, M, G, data, index_pat, coords, posterior_proba,
    Var_reg, family, offset_reg,
    formula_group, Var_group,
    prior_theta, prior_sigma, prior_proba,
    test.GR, seed,
    test.ICM, test.ICMregional, Wdist_LR, Wweight_SR, G_ICM, rho_ICM,
    export.predicteur)
}
\keyword{internal}
