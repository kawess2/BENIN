// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> 
#include <iostream> 
#include <Rmath.h> 
#include "Utilities_W.h"

using namespace Rcpp ;
using namespace std ;
using namespace arma ;

//// Fonctions de conversion vers R
// fct 1 : List calcGroupsCoords_cpp(const arma::mat& coords_NNA, const IntegerVector& index_NNA, const arma::mat& Neighborhood, IntegerVector& coords_max, 
//                          int& max_groups, bool& verbose)
// fct 2 : List calcGroupsW_cpp(const IntegerVector& W_i, const IntegerVector& W_p, const IntegerVector& subset, int& max_groups)
// fct 3 : List calcRadius_cpp(const arma::mat& coords, const NumericVector& sample, double& threshold, const LogicalVector& subset_bary, bool& verbose)
// fct 4 : List calcBlockW_cpp(const IntegerVector& W_i, const IntegerVector& W_p, const IntegerVector& site_order, 
//                    const NumericVector& dist_center, double dist_max, bool& verbose)

//  fct 1 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List calcGroupsCoords_cpp(const arma::mat& coords_NNA, const IntegerVector& index_NNA, int min_index_NNA, int max_index_NNA,
                          const arma::mat& Neighborhood, const IntegerVector& coords_max, 
                          int& max_groups, bool& verbose){
  
  List res = calcGroupsCoords_hpp(coords_NNA, 
                                  Rcpp::as < std::vector < int > >(index_NNA), 
                                  min_index_NNA,
                                  max_index_NNA,
                                  Neighborhood, 
                                  Rcpp::as < std::vector < int > >(coords_max), 
                                  max_groups, 
                                  verbose);
						  
  //  export
  return(res);  
}

// fct 2 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List calcGroupsW_cpp(const IntegerVector& W_i, const IntegerVector& W_p, const IntegerVector& subset, int& max_groups, bool verbose){
  //vector < int > W_i_cpp = W_i;
  //vector < int > W_p_cpp = W_p;
  
  List res = calcGroupsW_hpp(Rcpp::as < std::vector < int > >(W_i), 
                             Rcpp::as < std::vector < int > >(W_p), 
                             Rcpp::as < std::vector < int > >(subset), 
                             max_groups,
                             verbose);

    //  export
  return(res);
}

// fct 3 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List calcRadius_cpp(const arma::mat& coords, const NumericVector& sample, double& threshold, const LogicalVector& subset_bary, bool& verbose){
  
  List res = calcRadius_hpp(coords, 
                            Rcpp::as < std::vector < double > >(sample), 
                            threshold, 
                            Rcpp::as < std::vector < bool > >(subset_bary), 
                            verbose);
  
  // export
   return(res);     
}

// fct 4 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List calcBlockW_cpp(const IntegerVector& W_i, const IntegerVector& W_p, const IntegerVector& site_order, 
                    const NumericVector& dist_center, double dist_max, bool& verbose){

  List res = calcBlockW_hpp(Rcpp::as < std::vector < int > >(W_i), 
                            Rcpp::as < std::vector < int > >(W_p), 
                            Rcpp::as < std::vector < int > >(site_order), 
                            Rcpp::as < std::vector < double > >(dist_center), 
                            dist_max, 
                            verbose);
    
// export
return(res);     
}

  
