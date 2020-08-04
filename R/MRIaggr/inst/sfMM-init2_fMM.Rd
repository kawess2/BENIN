\name{init2_fMM}
\title{initialisation of the \code{\link{fMMalgo}} algorithm}
\alias{init2_fMM}
\description{
Initialisation of the \code{\link{fMMalgo}} algorithm : initialisation the distribution parameters for each group and of the regularization parameters.
For internal use.
}
\usage{
init2_fMM(n, M, G, data, Id,
    formula_reg, Var_reg, offset_reg,family,
    formula_group, Var_group, intercept_group,
    prior_theta, prior_sigma, prior_proba, posterior_proba,
    test.GR, seed,
    test.ICM, Wweight_SR, prior_prevalenceICM, rho_max, rho_ICM,
    test.ICMregional, distance_ref, Wdist_LR, coords, threshold, nbGroup_min, multiV,
    iter_max, trace)
}
\keyword{internal}
