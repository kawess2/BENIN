// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> 
#include <iostream> 
#include <Rmath.h> 
#include "Utilities_Potential.h"

using namespace Rcpp ;
using namespace std ;
using namespace arma ;

// fct 1 : List calcMultiPotential_cpp(const S4& W_SR, const S4& W_LR, const NumericVector& sample, double& threshold, const arma::mat& coords, 
//                                     NumericVector& distance_ref, int& nbGroup_min, bool& multiV, 
//                                     double& neutre)

// fct 1 : ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List calcMultiPotential_cpp(const S4& W_SR, const S4& W_LR, const NumericVector& sample, double& threshold, const arma::mat& coords, 
NumericVector& distance_ref, int& nbGroup_min, bool& multiV, 
double& neutre){ 
  
  List res = calcMultiPotential_hpp(W_SR, W_LR, 
                                    Rcpp::as < std::vector < double > >(sample), 
                                    threshold, coords, 
                                    Rcpp::as < std::vector < double > >(distance_ref), 
                                    nbGroup_min, multiV, neutre);

  return res;
}
