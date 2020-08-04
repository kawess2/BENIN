//#include <iostream>
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <iostream>
#include <progress.hpp>
#include "Utilities_General.hpp"
#include "Utilities_Potential.hpp"
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::depends("RcppProgress")]]

using namespace Rcpp ;
using namespace std ;
using namespace arma ;

//////////// Summary //////////////////////////////////////////
//> Main function (l 37)
// List calcPotts_cpp (const S4& Wdist_LR, const S4& Wweight_SR, arma::mat sample, NumericVector rho, const arma::mat& coords, bool type_reg,
//bool test_ICMregional, NumericVector distance_ref, double threshold, double neutre, int nbGroup_min, bool multiV, int type, bool distcat,
//bool prior, int iter_max, double critere, bool trace_iter, bool trace_radius)

//> Simulation functions (l 650)
// List simulPotts_cpp (const S4& Wweight_SR, const S4& Wdist_LR, arma::mat sample, const arma::mat& coords, NumericVector rho, NumericVector distance_ref, bool distcat, int iter_max, double critere, bool regional, bool trace)
// arma::mat simulPottsFast_cpp(const S4& Wweight_SR, arma::mat sample,  double rho, int n, int p, int iter_max)
////////////////////////////////////////////////////////


////> Main function ////
// [[Rcpp::export]]
List calcPotts_cpp (const S4& Wweight_SR, const S4& Wdist_LR, arma::mat sample, NumericVector rho, const arma::mat& coords, bool type_reg,
                    bool test_ICMregional, NumericVector distance_ref, double threshold, double neutre, int nbGroup_min, bool multiV, int type, bool distcat,
                    bool prior, int iter_max, double critere, bool trace_iter, bool trace_radius){
  
  // preparation
  int iter=1 ;
  int p=sample.n_cols;
  int n=sample.n_rows;
  arma::mat Wpred(n,p);
  arma::mat pred_global=sample;
  arma::mat pred_global_hist=sample;
  IntegerVector rang;
  int index_px;
  double norm;
  double val_critere=critere+1;
  NumericVector V(n);
  vector<double> sampleN(n);
  NumericVector Vsauve=V;
  double test_V;
  if(test_ICMregional==true){
    test_V=critere+1;
  }else{
    test_V=0; 
  }
  
  IntegerVector Wweight_SRi=Wweight_SR.slot("i");
  IntegerVector Wweight_SRp=Wweight_SR.slot("p");
  NumericVector Wweight_SRx=Wweight_SR.slot("x");
  List res_multipotentiel ;
  NumericVector rayon;
  
  if(prior==false){
    std::fill(sample.begin(),sample.end(),1);
  }
  
  
  while(iter <= iter_max && (val_critere > critere || test_V>critere)){
    
    rang = rank_hpp(runif(n)); // tirage de l ordre des sites
    
    // evaluation du potentiel regional
    if(test_ICMregional==true && (iter==1 || val_critere <= critere)){ 
      for(int iter_px=0 ; iter_px<n ; iter_px++) // recuperation du groupe necrose 
      {sampleN[iter_px] = pred_global(iter_px,p-1);}
      
      Vsauve = V;
      
      res_multipotentiel =  calcMultiPotential_hpp(Wweight_SR, Wdist_LR,  sampleN, threshold, coords, distance_ref, 
                                                   nbGroup_min, multiV, type,
                                                   neutre, trace_radius, distcat);    
      
      V = res_multipotentiel[0];
      //       if(type_reg){V = V*(V - sampleN);} // A CORRIGER
      //       test_V = max(abs(V-Vsauve)); // A CORRIGER
    }
    
    for(int iter_px=0 ; iter_px<n ; iter_px++){ // pour chaque pixel
      norm = 0.0;
      index_px = rang(iter_px)-1;
      
      for(int iter_p=0 ; iter_p<p ; iter_p++){ // pour chaque groupe
        Wpred(index_px,iter_p)=0.0;
        
        for(int iter_vois=Wweight_SRp[index_px]; iter_vois<Wweight_SRp[index_px+1]; iter_vois++){ // pour chaque voisin
          if(type_reg){
            Wpred(index_px,iter_p) = Wpred(index_px,iter_p) + Wweight_SRx[iter_vois]*pred_global(Wweight_SRi[iter_vois],iter_p)*max(0.0,pred_global(Wweight_SRi[iter_vois],iter_p)-sample(index_px,iter_p));
          }else{
            Wpred(index_px,iter_p) = Wpred(index_px,iter_p) + Wweight_SRx[iter_vois]*pred_global(Wweight_SRi[iter_vois],iter_p);
          }
        }
        
        // exponentielle rho
        if(test_ICMregional==false)
          {Wpred(index_px,iter_p) = exp(rho(0)*Wpred(index_px,iter_p)) ;}else{
            if(iter_p==p-1){
              Wpred(index_px,iter_p) = exp(rho(0)*Wpred(index_px,iter_p) + rho(1)*V(index_px));
            }else{
              Wpred(index_px,iter_p) = exp(rho(0)*Wpred(index_px,iter_p) + rho(1)*(1-V(index_px)));
            }
          }
          
          // calcul proba
          pred_global(index_px,iter_p) = sample(index_px,iter_p)*Wpred(index_px,iter_p);         
          
          norm = norm + pred_global(index_px,iter_p);
      } 
      
      // normalisation 
      for(int iter_p=0 ; iter_p<p ; iter_p++){
        pred_global(index_px,iter_p) = pred_global(index_px,iter_p)/norm;
      }
      
    }         
    
    val_critere=0;
    for(int iter_p=0 ; iter_p<p ; iter_p++){
      val_critere = max(val_critere,max(abs(pred_global_hist.col(iter_p)-pred_global.col(iter_p)))); // a passer en vectoriel
    }
    
    if(trace_iter==true){
      Rcout << "critere " << iter << " = " << val_critere << endl;
    }
    
    iter++;
    pred_global_hist = pred_global; 
  }  
  
  if(test_ICMregional==true){rayon = res_multipotentiel[1];}else{rayon = 0;}
  
  return Rcpp::List::create(Rcpp::Named("predit") = pred_global,
                            Rcpp::Named("pred_spatial") = Wpred,
                            Rcpp::Named("cv") = val_critere <= critere,
                                                           Rcpp::Named("rayon") = rayon,
                                                           Rcpp::Named("iter_max") = iter,
                                                           Rcpp::Named("V") = V);
  
}

////> Simulation functions ////
// [[Rcpp::export]]
List simulPotts_cpp (const S4& W_SR, const S4& W_LR, arma::mat sample, const arma::mat& coords,  
                     NumericVector rho, NumericVector distance_ref, bool distcat, 
                     int iter_max, double cv_criterion, bool regional, bool trace){
  
   // attention W_SR est lu par lignes donc si elle doit etre normalisee c est par lignes !!!
  
  //// initialization
  
  // diagnostic progress
  Progress testUser(trace*CST_PROGRESS, trace);
  double value_trace = 0 ;
  
  // variables
  int n=sample.n_rows;
  int p=sample.n_cols;
  vector<int> W_i=W_SR.slot("i");
  vector<int> W_p=W_SR.slot("p");
  vector<double> W_x=W_SR.slot("x");
  
  IntegerVector rang(n);
  int tirage_multinom[p];
  
  int index_px;
  arma::mat Wpred(n,p);
  double norm, proba_site[p];
  
  // regional
  arma::mat V(n,p);
  std::fill(V.begin(),V.end(),0);
  vector<double> sampleCol(n);
  List res_multipotentiel ;
  
  // diagnostic convergence
  bool check_cv=(cv_criterion>0);
  arma::mat proba_hist(check_cv*n,check_cv*p);
  double val_critere = cv_criterion + 1 ;
  bool test;
  
  //// main loop
  for(int iter=0; iter<iter_max ; iter++){
    
    // diagnostic
    if(trace && iter>=value_trace){
      value_trace = min(1.0*iter_max,value_trace + iter_max/CST_PROGRESS);
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
    rang = rank_hpp(runif(n));
    if(regional){   // regional potential
      for(int iter_p=0 ; iter_p<p ; iter_p++){
        
        for(int iter_obs=0 ; iter_obs<n ; iter_obs++){
          sampleCol[iter_obs]=sample(iter_obs,iter_p); // colvec to vector<double>
        }
        
        for(int iter_px=0 ; iter_px<n ; iter_px++){
          
          res_multipotentiel =  calcMultiPotential_hpp(W_SR, W_LR,  sampleCol, 0.1, 
                                                       coords, distance_ref, 
                                                       0, false, 2,
                                                       0.5, false, distcat);  
          
          V.col(iter_p) = as<arma::vec>(res_multipotentiel[0]);
        }
      }
    }

    for(int iter_px=0 ; iter_px<n ; iter_px++){ // pour chaque pixel
      
      norm = 0.0;
      index_px = rang[iter_px]-1;
      
      for(int iter_p=0 ; iter_p<p ; iter_p++){ // pour chaque groupe
        
        Wpred(index_px,iter_p)=0.0;
        
        // contribution de chaque voisin a la densite
        for(int iter_vois=W_p[index_px]; iter_vois<W_p[index_px+1]; iter_vois++){
          
          Wpred(index_px,iter_p) += W_x[iter_vois]*sample(W_i[iter_vois],iter_p);
          
        }
        
        // exponentielle rho
        Wpred(index_px,iter_p) = exp(rho[0]*Wpred(index_px,iter_p)+rho[1]*V(index_px,iter_p)) ;
        norm = norm + Wpred(index_px,iter_p);
      } 
      
      for(int iter_p=0; iter_p<p; iter_p++){
        proba_site[iter_p] = Wpred(index_px,iter_p)/norm;
        
        if(check_cv){
          if(iter>0){
            test = abs(proba_hist(iter_px,iter_p)-proba_site[iter_p]) ;
            if(test>val_critere){val_critere = test;}
            proba_hist(iter_px,iter_p) = proba_site[iter_p];  
          }
        }
      } 
      rmultinom(1, proba_site, p, tirage_multinom);
      
      for(int iter_p=0; iter_p<p; iter_p++){
        sample(index_px,iter_p) = tirage_multinom[iter_p];
      }
      
    }
    
    if(check_cv){
      if(val_critere<cv_criterion){break;}
    }
  }
  
    // cv
    bool test_cv;
    if(check_cv){
      test_cv = val_critere<cv_criterion;
    }else{
      test_cv =NA_REAL;
    }
    
    // export
    return Rcpp::List::create(Rcpp::Named("simulation") = sample,
                              Rcpp::Named("V") = V,
                              Rcpp::Named("cv_criterion") = cv_criterion,
                              Rcpp::Named("cv") = test_cv
    );
}

// // [[Rcpp::export]]
// arma::mat simulPottsFast_cpp(const IntegerVector& W_i, const IntegerVector& W_p, const NumericVector& W_x, 
//                              arma::mat sample,  double rho, int n, int p,int iter_nb){
//  
//   // attention W est lu par lignes donc si elle doit etre normalisee c est par lignes !!!
//   int index_px;
//   double norm;   
//   arma::mat Wpred(n,p);
//   IntegerVector rang(n);
//   int tirage_multinom[p];
//   double proba_site[p];
//   
//   for(int iter=0 ; iter<iter_nb ; iter++){
//     
//     rang = rank_hpp(runif(n)); // tirer aleatoirement l ordre des sites
//  
//     for(int iter_px=0 ; iter_px<n ; iter_px++){ // pour chaque pixel
//       
//       norm = 0.0;
//       index_px = rang(iter_px)-1;
//       
//       for(int iter_p=0 ; iter_p<p ; iter_p++){ // pour chaque groupe
//         
//         Wpred(index_px,iter_p)=0.0;
//         
//         // contribution de chaque voisin a la densite
//         for(int iter_vois=W_p[index_px]; iter_vois<W_p[index_px+1]; iter_vois++){
//           Wpred(index_px,iter_p) = Wpred(index_px,iter_p) + W_x[iter_vois]*sample(W_i[iter_vois],iter_p);
//         }
//         
//         // exponentielle rho
//         Wpred(index_px,iter_p) = exp(rho*Wpred(index_px,iter_p)) ;
//         norm = norm + Wpred(index_px,iter_p);
//       } 
//       
//       for(int iter_p=0; iter_p<p; iter_p++){
//         proba_site[iter_p] = Wpred(index_px,iter_p)/norm;
//       } 
//       
//       rmultinom(1, proba_site, p, tirage_multinom);
//       
//       for(int iter_p=0; iter_p<p; iter_p++){
//         sample(index_px,iter_p) = tirage_multinom[iter_p];
//       }
//       
//     }
//   }
//   
//   
//   return(sample);
// }
