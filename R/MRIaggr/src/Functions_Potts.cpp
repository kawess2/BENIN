// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]
#include <RcppArmadillo.h> 
#include <progress.hpp> 
#include <Rmath.h> 
#include <iostream> 
#include "Utilities_General.h"
#include "Utilities_Potential.h"

using namespace Rcpp ;
using namespace std ;
using namespace arma ;

//////////// Summary //////////////////////////////////////////
//> Main function
// List calcPotts_cpp (const S4& W_SR, const S4& W_LR, arma::mat sample, NumericVector rho, const arma::mat& coords, 
//                    const IntegerVector& site_order, int iter_max, double cv_criterion, 
//                    bool test_regional, NumericVector distance_ref, double threshold, double neutre, int nbGroup_min, bool multiV, bool last_vs_others, 
//                    bool prior, bool type_reg, int verbose)

//> Simulation functions
// List simulPotts_cpp (const S4& W_SR, const S4& W_LR, arma::mat sample, const arma::mat& coords, 
//                     const IntegerVector& site_order, 
//                     NumericVector rho, NumericVector distance_ref, 
//                     int iter_max, double cv_criterion, bool regional, bool verbose)
// arma::mat simulPottsFast_cpp(const IntegerVector& W_i, const IntegerVector& W_p, const NumericVector& W_x, 
//                             const IntegerVector& site_order, 
//                             arma::mat sample, double rho, int n, int p, int iter_nb)
////////////////////////////////////////////////////////

////> Main function ////
// [[Rcpp::export]]
List calcPotts_cpp (const S4& W_SR, const S4& W_LR, arma::mat sample, NumericVector rho, const arma::mat& coords, 
                    const IntegerVector& site_order, int iter_max, double cv_criterion, 
                    bool test_regional, NumericVector distance_ref, double threshold, double neutre, int nbGroup_min, bool multiV, bool last_vs_others, 
                    bool prior, bool type_reg, int verbose){
  
  //// preparation
  int p = sample.n_cols;
  int n = sample.n_rows;
  
  IntegerVector W_SRi = W_SR.slot("i");
  IntegerVector W_SRp = W_SR.slot("p");
  NumericVector W_SRx = W_SR.slot("x");
  
  arma::mat Wpred(n, p);
   
  vector < vector < double > > pred_global(p);
  for(int iter_p = 0 ; iter_p < p ; iter_p ++){
    pred_global[iter_p].resize(n);
    for(int iter_px = 0 ; iter_px < n ; iter_px++){
      pred_global[iter_p][iter_px] = sample(iter_px, iter_p);
    }
  }
  vector < vector < double > > pred_global_hist = pred_global;
    
  if(prior == false){
    std::fill(sample.begin(), sample.end(), 1);
  }
  
  arma::mat V(n, p);
   
  IntegerVector rang;
  bool no_site_order = (site_order[0] < 0);
  if(no_site_order == false){
    rang = site_order;
  }
  
  int index_px;
  double norm;
  double val_criterion = cv_criterion + 1;
  double diff_criterion = 0, diff; 
  double cv_criterion2 = n * cv_criterion ; 
  int iter_updateV = 0 ;
  
  vector < double > res_multipotentiel(n) ;
  int iter = 0 ;
  
  //// estimation
  
  while(iter < iter_max && ( (val_criterion > cv_criterion) || (test_regional == true && iter_updateV != iter) ) ){

    iter++;
    pred_global_hist = pred_global; 
    
    if(no_site_order){
      rang = rank_hpp(runif(n)) - 1; // tirer aleatoirement l ordre des sites
    }
   
    //// regional potential
    if(test_regional == true){
      if(iter == 1 || val_criterion <= cv_criterion || diff_criterion > cv_criterion2){ 
      iter_updateV = iter ;
      
      for(int iter_p = (p - 1) ; iter_p>= 0 ; iter_p --){

        if(last_vs_others == false || iter_p == (p - 1)){
        
        res_multipotentiel =  calcMultiPotential_hpp(W_SR, W_LR, pred_global[iter_p], 
                                                     threshold, coords, Rcpp::as < std::vector < double > >(distance_ref), 
                                                     nbGroup_min, multiV, 
                                                     neutre)[0]; 
            
          for(int iter_px = 0 ; iter_px < n ; iter_px++){
           
            V(iter_px, iter_p) = res_multipotentiel[iter_px];
            
            if(type_reg){
              V(iter_px, iter_p) = V(iter_px, iter_p)  * (V(iter_px, iter_p) - pred_global[iter_p][iter_px]);
            }
          }
        
        }else{
          V.col(iter_p) = 1 - V.col(p - 1);
        }
      }
      
      }else{
        cv_criterion2 /= 1.1;
      } 
    }

    //// update site probabilities
    val_criterion = 0;
    diff_criterion = 0;
    
    for(int iter_px = 0 ; iter_px < n ; iter_px++){ // pour chaque pixel
      
      norm = 0.0;
      index_px = rang(iter_px);
      
      for(int iter_p = 0 ; iter_p < p ; iter_p++){ // pour chaque groupe
        Wpred(index_px, iter_p) = 0.0;
        
        for(int iter_vois = W_SRp[index_px]; iter_vois < W_SRp[index_px + 1]; iter_vois++){ // pour chaque voisin
          if(type_reg){
            Wpred(index_px, iter_p) += W_SRx[iter_vois]*pred_global[iter_p][W_SRi[iter_vois]]*max(0.0, pred_global[iter_p][W_SRi[iter_vois]]-sample(index_px, iter_p));
            }else{
            Wpred(index_px, iter_p) += W_SRx[iter_vois] *pred_global[iter_p][W_SRi[iter_vois]];
          }
        }
        
        // exponentielle rho
        if(test_regional == false){
          pred_global[iter_p][index_px] = sample(index_px, iter_p) *exp(rho(0) *Wpred(index_px, iter_p)) ;
        }else{
          pred_global[iter_p][index_px] = sample(index_px, iter_p) *exp(rho(0) *Wpred(index_px, iter_p) + rho(1) *V(index_px, iter_p));
        } 
        norm += pred_global[iter_p][index_px];
      } 
      
      // normalisation 
      for(int iter_p = 0 ; iter_p < p ; iter_p++){
        pred_global[iter_p][index_px] /= norm;
        diff = abs(pred_global_hist[iter_p][index_px]-pred_global[iter_p][index_px]);
        diff_criterion += diff;
        val_criterion = max(val_criterion, diff);
      }
      
    }         

    if(verbose > 0){
      if(verbose == 1){Rcout << "*" ;}
      if(verbose == 2){Rcout << "iteration " << iter << " : totaldiff = " << diff_criterion << " | maxdiff = " << val_criterion << endl;}
    }

  }  
  
  if(verbose == 1){Rcout << endl;}
  
  //// export
  arma::mat spatialPrior(n, p);
  for(int iter_px = 0 ; iter_px < n ; iter_px++){
    index_px = rang(iter_px);
    for(int iter_p = 0 ; iter_p < p ; iter_p++){
      if(test_regional == false){
        spatialPrior(index_px, iter_p) = exp(rho(0) *Wpred(index_px, iter_p)) ;
      }else{
          spatialPrior(index_px, iter_p) = exp(rho(0) *Wpred(index_px, iter_p) + rho(1) *V(index_px, iter_p));
      } 
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("predicted") = pred_global, 
                            Rcpp::Named("spatialPrior") = spatialPrior, 
                            Rcpp::Named("cv") = val_criterion <= cv_criterion, 
                            Rcpp::Named("iter_max") = iter
                            );
  
}

////> Simulation functions ////
// [[Rcpp::export]]
List simulPotts_cpp (const S4& W_SR, const S4& W_LR, arma::mat sample, const arma::mat& coords, 
                     const IntegerVector& site_order, 
                     NumericVector rho, NumericVector distance_ref, 
                     int iter_max, double cv_criterion, bool regional, bool verbose){
  
   // attention W_SR est lu par lignes donc si elle doit etre normalisee c est par lignes !!!
  
  //// initialization
  
  // diagnostic progress
  Progress testUser(verbose *CST_PROGRESS, verbose);
  double value_trace = 0 ;
  
  // variables
  int n = sample.n_rows;
  const int p = sample.n_cols;
  vector < int > W_i = W_SR.slot("i");
  vector < int > W_p = W_SR.slot("p");
  vector < double > W_x = W_SR.slot("x");
  
  IntegerVector rang(n);
  bool no_site_order = (site_order[0] < 0);
  if(no_site_order == false){
    rang = site_order;
  }
  IntegerVector tirage_multinom(p); // int tirage_multinom[p]; //  (generate a warning on linux : warning variable length arrays are a C99 feature)
  NumericVector proba_site(p); //double proba_site[p]; // (generate a warning on linux : warning variable length arrays are a C99 feature)
  
  int index_px;
  arma::mat Wpred(n, p);
  double norm;
  
  // regional
  arma::mat V(n, p);
  std::fill(V.begin(), V.end(), 0);
  vector < double > sampleCol(n);
  List res_multipotentiel ;
  
  // diagnostic convergence
  bool check_cv = (cv_criterion > 0);
  arma::mat proba_hist(check_cv *n, check_cv *p);
  double val_criterion = cv_criterion + 1 ;
  bool test;
  
  //// main loop
  for(int iter = 0; iter < iter_max ; iter++){
    
    // diagnostic
    if(verbose && iter>=value_trace){
      value_trace = min(1.0 *iter_max, value_trace + iter_max / CST_PROGRESS);
      testUser.increment();
    }
    if (Progress::check_abort() ){
      sample.fill(NA_REAL);
      V.fill(NA_REAL);
      return Rcpp::List::create(Rcpp::Named("simulation") = sample, 
                                Rcpp::Named("V") = V, 
                                Rcpp::Named("cv") = false
      );
    }

    // sort site order
    if(no_site_order){
      rang = rank_hpp(runif(n)) - 1; // tirer aleatoirement l ordre des sites
    }
    
    if(regional){   // regional potential
      for(int iter_p = 0 ; iter_p < p ; iter_p++){
        
        for(int iter_obs = 0 ; iter_obs < n ; iter_obs++){
          sampleCol[iter_obs] = sample(iter_obs, iter_p); // colvec to vector < double > 
        }
        
        for(int iter_px = 0 ; iter_px < n ; iter_px++){
          
          res_multipotentiel =  calcMultiPotential_hpp(W_SR, W_LR,  sampleCol, 0.01, 
                                                       coords, Rcpp::as < std::vector < double > >(distance_ref), 
                                                       true, 10, 0.5);  
          
          V.col(iter_p) = as < arma::vec >(res_multipotentiel[0]);
        }
      }
    }

    for(int iter_px = 0 ; iter_px < n ; iter_px++){ // pour chaque pixel
      
      norm = 0.0;
      index_px = rang[iter_px];
      
      for(int iter_p = 0 ; iter_p < p ; iter_p++){ // pour chaque groupe
        
        Wpred(index_px, iter_p) = 0.0;
        
        // contribution de chaque voisin a la densite
        for(int iter_vois = W_p[index_px]; iter_vois < W_p[index_px + 1]; iter_vois++){
          
          Wpred(index_px, iter_p) += W_x[iter_vois] *sample(W_i[iter_vois], iter_p);
          
        }
        
        // exponentielle rho
        Wpred(index_px, iter_p) = exp(rho[0] * Wpred(index_px, iter_p) + rho[1] *V(index_px, iter_p)) ;
        norm = norm + Wpred(index_px, iter_p);
      } 
      
      for(int iter_p = 0; iter_p < p; iter_p++){
        proba_site[iter_p] = Wpred(index_px, iter_p) / norm;
        
        if(check_cv){
          if(iter > 0){
            test = abs(proba_hist(iter_px, iter_p) - proba_site[iter_p]) ;
            if(test > val_criterion){val_criterion = test;}
            proba_hist(iter_px, iter_p) = proba_site[iter_p];  
          }
        }
      } 
      rmultinom(1, proba_site.begin(), p, tirage_multinom.begin()); // rmultinom(1, proba_site, p, tirage_multinom); //  (alternative version with the warning)
      
      for(int iter_p = 0; iter_p < p; iter_p++){
        sample(index_px, iter_p) = tirage_multinom[iter_p];
      }
      
    }
    
    if(check_cv){
      if(val_criterion < cv_criterion){break;}
    }
  }
  
    // cv
    bool test_cv;
    if(check_cv){
      test_cv = val_criterion < cv_criterion;
    }else{
      test_cv  = NA_REAL;
    }
    
    // export
    return Rcpp::List::create(Rcpp::Named("simulation") = sample, 
                              Rcpp::Named("V") = V, 
                              Rcpp::Named("cv_criterion") = cv_criterion, 
                              Rcpp::Named("cv") = test_cv
    );
}

// [[Rcpp::export]]
arma::mat simulPottsFast_cpp(const IntegerVector& W_i, const IntegerVector& W_p, const NumericVector& W_x, 
                             const IntegerVector& site_order, 
                             arma::mat sample, double rho, int n, const int p, int iter_nb){
 
  // attention W est lu par lignes donc si elle doit etre normalisee c est par lignes !!!
  int index_px;
  double norm;
  arma::mat Wpred(n, p);
  IntegerVector rang(n);

  IntegerVector tirage_multinom(p); 
  NumericVector proba_site(p);
//  int tirage_multinom[p]; //// (generate a warning on linux : warning variable length arrays are a C99 feature)
//  double proba_site[p]; // // (generate a warning on linux : warning variable length arrays are a C99 feature)
  
//// sinon essayer une allocation dynamic
//   proba_site=malloc(p*sieof(double));
//   int *tirage_multinom;
//   tirage_multinom=calloc(p,sizeof(int));
//   double *proba_site;
//   proba_site=calloc(p,sizeof(double));

  bool no_site_order = (site_order[0] < 0);
  if(no_site_order == false){
    rang = site_order;
  }

  for(int iter = 0 ; iter < iter_nb ; iter++){
    
    if(no_site_order){
      rang = rank_hpp(runif(n)) - 1; // tirer aleatoirement l ordre des sites
      }
 
    for(int iter_px = 0 ; iter_px < n ; iter_px++){ // pour chaque pixel
      
      norm = 0.0;
      index_px = rang(iter_px);
      
      for(int iter_p = 0 ; iter_p < p ; iter_p++){ // pour chaque groupe
        
        Wpred(index_px, iter_p) = 0.0;
        
        // contribution de chaque voisin a la densite
        for(int iter_vois = W_p[index_px]; iter_vois < W_p[index_px + 1]; iter_vois++){
          Wpred(index_px, iter_p) += W_x[iter_vois] *sample(W_i[iter_vois], iter_p);
        }
        
        // exponentielle rho
        Wpred(index_px, iter_p) = exp(rho *Wpred(index_px, iter_p)) ;
        norm = norm + Wpred(index_px, iter_p);
      } 
      
      for(int iter_p = 0; iter_p < p; iter_p++){
        proba_site[iter_p] = Wpred(index_px, iter_p) / norm;
      } 
      
      rmultinom(1, proba_site.begin(), p, tirage_multinom.begin());
      // rmultinom(1, proba_site, p, tirage_multinom);   //  (alternative version with the warning)
      
      for(int iter_p = 0; iter_p < p; iter_p++){
        sample(index_px, iter_p) = tirage_multinom[iter_p];
      }
      
    }
  }
  
  //// sinon essayer une allocation dynamic
//   free(proba_site);
//   free(tirage_multinom);
  
  
  return(sample);
}
