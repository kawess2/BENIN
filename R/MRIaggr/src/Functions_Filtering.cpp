// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> 
#include <iostream> 
#include <Rmath.h> 
#include "Utilities_General.h"

using namespace Rcpp ;
using namespace std ;
using namespace arma ;


// arma::mat filtrage2D_cpp(const arma::mat& M_data, const arma::mat& M_operateur, const arma::mat& index_data, 
//                    bool bilateral, bool na_rm)
// List filtrage3D_cpp(const NumericVector& Vec_data, const IntegerVector& p_data, 
//                    const NumericVector& Vec_operateur, const IntegerVector& p_operateur, 
//                    const arma::mat& index_data, bool bilateral, bool na_rm)
// arma::mat filtrage2Dmed_cpp(const arma::mat& M_data, const arma::mat& M_operateur, const arma::mat& index_data, 
//                            bool na_rm)
// arma::mat filtrage3Dmed_cpp(const NumericVector& Vec_data, const IntegerVector& p_data, 
//                       const NumericVector& Vec_operateur, const IntegerVector& p_operateur, 
//                       const arma::mat& index_data, bool na_rm)

//  1 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List filtrage2D_cpp(const arma::mat& M_data, const arma::mat& M_operateur, const arma::mat& index_data, 
                    bool bilateral, bool na_rm){

  int const n_data = index_data.n_rows;
  vector < int > p_data(2);
  p_data[0] = M_data.n_rows ; p_data[1] = M_data.n_cols ;
  vector < int > p_operateur(2);
  p_operateur[0] = M_operateur.n_rows ; p_operateur[1] = M_operateur.n_cols;

  vector < int > v0(p_operateur[0]);
  for(int iter = 0 ; iter < p_operateur[0] ; iter++){v0[iter] = iter;}
  vector < int > v1(p_operateur[1]);
  for(int iter = 0 ; iter < p_operateur[1] ; iter++){v1[iter] = iter;}
  vector < double > p_ref(2); // vector < int > p_ref(2);
  p_ref[0] = medianVecInteger_hpp(v0) ; p_ref[1] = medianVecInteger_hpp(v1); 

  NumericVector M_vector(n_data);
  for(int iter_px = 0 ; iter_px < n_data ; iter_px++)
  {M_vector(iter_px) = M_data(index_data(iter_px, 0), index_data(iter_px, 1)) ;}
  
  double const sigma = sqrt(var(M_vector));
  double const dnorm_sigma = R::dnorm(0.0, 0.0, sigma, false); 
 
 
  arma::mat Mres(p_data[0], p_data[1]);  
  std::fill(Mres.begin(), Mres.end(), NA_REAL);
  arma::mat Wres(p_data[0], p_data[1]);  
  std::fill(Wres.begin(), Wres.end(), NA_REAL);
  bool test1a, test1b, test2a, test2b, break_loop ; 
  double tempo, sumNNA, Mintensite, Mvoisin ;
  
  for(int iter_px = 0 ; iter_px < n_data ; iter_px++){
    // initialisation
    Mres(index_data(iter_px, 0), index_data(iter_px, 1)) = 0;
    sumNNA = 0;
    break_loop = false;
      
    // mise en place des matrices
    for(int iter_l = 0 ; iter_l < p_operateur[0] ; iter_l++){

      for(int iter_c = 0 ; iter_c < p_operateur[1] ; iter_c++){   

         Mintensite = 1;
         Mvoisin = NA_REAL;
        
         test1a = index_data(iter_px, 0) + iter_l-p_ref[0] >= 0;
         test1b = index_data(iter_px, 0) + iter_l-p_ref[0] < p_data[0];
         test2a = index_data(iter_px, 1) + iter_c-p_ref[1] >= 0;
         test2b = index_data(iter_px, 1) + iter_c-p_ref[1] < p_data[1];
    
    if(test1a * test1b * test2a * test2b == 1 && R_IsNA(M_data(index_data(iter_px, 0) + iter_l-p_ref[0], index_data(iter_px, 1) + iter_c-p_ref[1])) == 0){      
      // elements voisins
      if(R_IsNA(M_operateur(iter_l, iter_c)) == 0){
        Mvoisin = M_data(index_data(iter_px, 0) + iter_l-p_ref[0], 
        index_data(iter_px, 1) + iter_c-p_ref[1]);
        
        // intensite des pixels
        if(bilateral == 1 && M_operateur(iter_l, iter_c) != 0){
          tempo = Mvoisin-M_data(index_data(iter_px, 0), index_data(iter_px, 1));
          Mintensite = R::dnorm(tempo, 0.0, sigma, false) / dnorm_sigma;
        }
        
      }
    }
     
    // ajustement des NA et maj
    if(R_IsNA(Mvoisin) == 0){
      sumNNA += abs(M_operateur(iter_l, iter_c) * Mintensite);
 
      Mres(index_data(iter_px, 0), index_data(iter_px, 1)) += M_operateur(iter_l, iter_c) * Mintensite * Mvoisin;
    }else{
      
      if(na_rm){
        Mres(index_data(iter_px, 0), index_data(iter_px, 1)) = NA_REAL;
        break_loop = true;
        break;
      }
            
    }
    }
    if(break_loop){break;}
    }    

    // correction des NA
    Wres(index_data(iter_px, 0), index_data(iter_px, 1)) = sumNNA;
//    if(break_loop == false && sumNNA > 0)
//    {Mres(index_data(iter_px, 0), index_data(iter_px, 1)) = Mres(index_data(iter_px, 0), index_data(iter_px, 1)) / sumNNA;}
 }  

    return(List::create(Named("Mres")  = Mres, 
                        Named("Wres")  = Wres                          
                ));
}

//  2 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List filtrage3D_cpp(const NumericVector& Vec_data, const IntegerVector& p_data, 
                    const NumericVector& Vec_operateur, const IntegerVector& p_operateur, 
                    const arma::mat& index_data, bool bilateral, bool na_rm){


  arma::cube A_data(Vec_data.begin(), p_data[0], p_data[1], p_data[2]);
  arma::cube A_operateur(Vec_operateur.begin(), p_operateur[0], p_operateur[1], p_operateur[2]);
  int const n_data = index_data.n_rows;
  
  vector < int > v0(p_operateur[0]);
  for(int iter = 0 ; iter < p_operateur[0] ; iter++){v0[iter] = iter;}
  vector < int > v1(p_operateur[1]);
  for(int iter = 0 ; iter < p_operateur[1] ; iter++){v1[iter] = iter;}
  vector < int > v2(p_operateur[2]);
  for(int iter = 0 ; iter < p_operateur[2] ; iter++){v2[iter] = iter;}
  vector < double > p_ref(3); // vector < int > p_ref(3);
  p_ref[0] = medianVecInteger_hpp(v0) ; p_ref[1] = medianVecInteger_hpp(v1);  p_ref[2] = medianVecInteger_hpp(v2);

  // variance des donnees pour le lissage
  NumericVector M_vector(n_data);
  for(int iter_px = 0 ; iter_px < n_data ; iter_px++)
  {M_vector(iter_px) = A_data(index_data(iter_px, 0), index_data(iter_px, 1), index_data(iter_px, 2)) ;}
  
  double const sigma = sqrt(var(M_vector));
  double const dnorm_sigma = R::dnorm(0.0, 0.0, sigma, false); 
 
  // stockage resultats
  arma::cube Mres(p_data(0), p_data(1), p_data(2));
  std::fill(Mres.begin(), Mres.end(), NA_REAL);
  arma::cube Wres(p_data(0), p_data(1), p_data(2));
  std::fill(Wres.begin(), Wres.end(), NA_REAL);
  
  bool test1a, test1b, test2a, test2b, test3a, test3b, break_loop ; 
  double tempo, sumNNA, Mvoisin, Mintensite;
    
  for(int iter_px = 0 ; iter_px < n_data ; iter_px++){
  
    // initialisation
    Mres(index_data(iter_px, 0), index_data(iter_px, 1), index_data(iter_px, 2)) = 0;
    sumNNA = 0;
    break_loop = false;
       
    // mise en place des matrices
    for(int iter_l = 0 ; iter_l < p_operateur[0] ; iter_l++){
      
      for(int iter_c = 0 ; iter_c < p_operateur[1] ; iter_c++){      
        
         for(int iter_d = 0 ; iter_d < p_operateur[2] ; iter_d++){ 
           Mintensite = 1;
           Mvoisin = NA_REAL;
    
         test1a = index_data(iter_px, 0) + iter_l-p_ref[0] >= 0;
         test1b = index_data(iter_px, 0) + iter_l-p_ref[0] < p_data[0];
         test2a = index_data(iter_px, 1) + iter_c-p_ref[1] >= 0;
         test2b = index_data(iter_px, 1) + iter_c-p_ref[1] < p_data[1];
         test3a = index_data(iter_px, 2) + iter_d-p_ref[2] >= 0;
         test3b = index_data(iter_px, 2) + iter_d-p_ref[2] < p_data[2];
      
    
    if(test1a * test1b * test2a * test2b * test3a * test3b == 1 && R_IsNA(A_data(index_data(iter_px, 0) + iter_l-p_ref[0], index_data(iter_px, 1) + iter_c-p_ref[1], index_data(iter_px, 2) + iter_d-p_ref[2])) == 0){
      // elements voisins
       if(R_IsNA(A_operateur(iter_l, iter_c, iter_d)) == 0){
              Mvoisin = A_data(index_data(iter_px, 0) + iter_l-p_ref[0], 
                               index_data(iter_px, 1) + iter_c-p_ref[1], 
                               index_data(iter_px, 2) + iter_d-p_ref[2]);
            
      // intensite des pixels
      if(bilateral == 1){
      tempo = Mvoisin-A_data(index_data(iter_px, 0), index_data(iter_px, 1), index_data(iter_px, 2));
      Mintensite = R::dnorm(tempo, 0.0, sigma, false) / dnorm_sigma;
      }
       
      }
      
    }

    // ajustement des NA et maj
    if(R_IsNA(Mvoisin) == 0){
        sumNNA += abs(A_operateur(iter_l, iter_c, iter_d) * Mintensite);

        Mres(index_data(iter_px, 0), index_data(iter_px, 1), index_data(iter_px, 2)) += A_operateur(iter_l, iter_c, iter_d) * Mintensite * Mvoisin;
                     
    }else{
      
      if(na_rm)
      { Mres(index_data(iter_px, 0), index_data(iter_px, 1), index_data(iter_px, 2)) = NA_REAL;
        break_loop = true;
        break;
      }
            
    }
        }
         if(break_loop){break;}
    }
    if(break_loop){break;}
    }    

    // correction des NA
    Wres(index_data(iter_px, 0), index_data(iter_px, 1), index_data(iter_px, 2)) = sumNNA;
//    if(break_loop == false && sumNNA > 0)
//    {Mres(index_data(iter_px, 0), index_data(iter_px, 1), index_data(iter_px, 2)) = Mres(index_data(iter_px, 0), index_data(iter_px, 1), index_data(iter_px, 2)) / sumNNA;}

 }  

    return(List::create(Named("Mres")  = Mres, 
                        Named("Wres")  = Wres                          
                ));
}

//  3 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat filtrage2Dmed_cpp(const arma::mat& M_data, const arma::mat& M_operateur, const arma::mat& index_data, 
                            bool na_rm){
  
  int const n_data = index_data.n_rows;
  vector < int > p_data(2);
  p_data[0] = M_data.n_rows ; p_data[1] = M_data.n_cols ;
  vector < int > p_operateur(2);
  p_operateur[0] = M_operateur.n_rows ; p_operateur[1] = M_operateur.n_cols;

  vector < int > v0(p_operateur[0]);
  for(int iter = 0 ; iter < p_operateur[0] ; iter++){v0[iter] = iter;}
  vector < int > v1(p_operateur[1]);
  for(int iter = 0 ; iter < p_operateur[1] ; iter++){v1[iter] = iter;}
  vector < double > p_ref(2); //vector < int > p_ref(2);
  p_ref[0] = medianVecInteger_hpp(v0) ; p_ref[1] = medianVecInteger_hpp(v1);
  vector < double > Mvoisin ;

  arma::mat Mres(p_data[0], p_data[1]);
  std::fill(Mres.begin(), Mres.end(), NA_REAL);
  bool test1a, test1b, test2a, test2b, break_loop ; 
 
  for(int iter_px = 0 ; iter_px < n_data ; iter_px++){
    
    // initialisation
    Mres(index_data(iter_px, 0), index_data(iter_px, 1)) = 0;
    Mvoisin.resize(0);Mvoisin.reserve(p_ref[0]*p_ref[1]);
    break_loop = false;
      
      
    // mise en place des matrices
    for(int iter_l = 0 ; iter_l < p_operateur[0] ; iter_l++){
      
      for(int iter_c = 0 ; iter_c < p_operateur[1] ; iter_c++){      
         test1a = index_data(iter_px, 0) + iter_l-p_ref[0] >= 0;
         test1b = index_data(iter_px, 0) + iter_l-p_ref[0] < p_data[0];
         test2a = index_data(iter_px, 1) + iter_c-p_ref[1] >= 0;
         test2b = index_data(iter_px, 1) + iter_c-p_ref[1] < p_data[1];
                          
    // elements voisins
    if( test1a * test1b * test2a * test2b == 1 && R_IsNA(M_data(index_data(iter_px, 0) + iter_l-p_ref[0], index_data(iter_px, 1) + iter_c-p_ref[1])) == 0){
          
      if(R_IsNA(M_operateur(iter_l, iter_c)) == 0  && M_operateur(iter_l, iter_c) != 0){
      Mvoisin.push_back(M_data(index_data(iter_px, 0) + iter_l-p_ref[0], 
                          index_data(iter_px, 1) + iter_c-p_ref[1]));
      }
      
    }else{
      if(na_rm){
        Mres(index_data(iter_px, 0), index_data(iter_px, 1))= NA_REAL;
        break_loop = true;
        break;                       
      }
    }
        
    }
    if(break_loop){break;}
    }    
  
    // calcul de la mediane
    if(R_IsNA(Mres(index_data(iter_px, 0), index_data(iter_px, 1))) == 0){
    Mres(index_data(iter_px, 0), index_data(iter_px, 1)) = medianVecDouble_hpp(Mvoisin);
    }
    
  
 }  

 return Mres;
}

//  4 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
NumericVector filtrage3Dmed_cpp(const NumericVector& Vec_data, const IntegerVector& p_data, 
                       const NumericVector& Vec_operateur, const IntegerVector& p_operateur, 
                       const arma::mat& index_data, bool na_rm){
  
    // conversion en array
   arma::cube A_data(Vec_data.begin(), p_data[0], p_data[1], p_data[2]);
   arma::cube A_operateur(Vec_operateur.begin(), p_operateur[0], p_operateur[1], p_operateur[2]);
   int const n_data = index_data.n_rows;

  vector < int > v0(p_operateur[0]);
  for(int iter = 0 ; iter < p_operateur[0] ; iter++){v0[iter] = iter;}
  vector < int > v1(p_operateur[1]);
  for(int iter = 0 ; iter < p_operateur[1] ; iter++){v1[iter] = iter;}
  vector < int > v2(p_operateur[2]);
  for(int iter = 0 ; iter < p_operateur[2] ; iter++){v2[iter] = iter;}

  vector < double > p_ref(3);  //vector < int > p_ref(3);
   p_ref[0] = medianVecInteger_hpp(v0) ; p_ref[1] = medianVecInteger_hpp(v1);  p_ref[2] = medianVecInteger_hpp(v2);

  vector < double > Mvoisin ;
 
  arma::cube Mres(p_data[0], p_data[1], p_data[2]);
  std::fill(Mres.begin(), Mres.end(), NA_REAL);
    
    bool test1a, test1b, test2a, test2b, test3a, test3b, break_loop ; 
   
      for(int iter_px = 0 ; iter_px < n_data ; iter_px++){
          
            // initialisation
          Mres(index_data(iter_px, 0), index_data(iter_px, 1), index_data(iter_px, 2)) = 0;
          Mvoisin.resize(0);Mvoisin.reserve(p_ref[0]*p_ref[1]*p_data[2]);
          break_loop = false;
          
            // mise en place des matrices
            for(int iter_l = 0 ; iter_l < p_operateur[0] ; iter_l++){
               
                for(int iter_c = 0 ; iter_c < p_operateur[1] ; iter_c++){      
                  
                    for(int iter_d = 0 ; iter_d < p_operateur[2] ; iter_d++){ 
                       
                     test1a = index_data(iter_px, 0) + iter_l-p_ref[0] >= 0;
                     test1b = index_data(iter_px, 0) + iter_l-p_ref[0] < p_data[0];
                     test2a = index_data(iter_px, 1) + iter_c-p_ref[1] >= 0;
                     test2b = index_data(iter_px, 1) + iter_c-p_ref[1] < p_data[1];
                     test3a = index_data(iter_px, 2) + iter_d-p_ref[2] >= 0;
                     test3b = index_data(iter_px, 2) + iter_d-p_ref[2] < p_data[2];
                                      
                  // elements voisins
                if( test1a * test1b * test2a * test2b * test3a * test3b == 1 && R_IsNA(A_data(index_data(iter_px, 0) + iter_l-p_ref[0], index_data(iter_px, 1) + iter_c-p_ref[1], index_data(iter_px, 2) + iter_d-p_ref[2])) == 0){
                        
                      if(R_IsNA(A_operateur(iter_l, iter_c, iter_d)) == 0  && A_operateur(iter_l, iter_c, iter_d) != 0){
                        Mvoisin.push_back(A_data(index_data(iter_px, 0) + iter_l-p_ref[0], 
                                            index_data(iter_px, 1) + iter_c-p_ref[1], 
                                            index_data(iter_px, 2) + iter_d-p_ref[2]));
                        }
                    
                    }else{
                        if(na_rm){
                            Mres(index_data(iter_px, 0), index_data(iter_px, 1), index_data(iter_px, 2))= NA_REAL;
                            break_loop = true;
                            break;                       
                          }
                      }
                     }                    
                    if(break_loop){break;}
                   }                  
            if(break_loop){break;}
            }   
            
        
            // calcul de la mediane
          if(R_IsNA(Mres(index_data(iter_px, 0), index_data(iter_px, 1), index_data(iter_px, 2))) == 0)
            {Mres(index_data(iter_px, 0), index_data(iter_px, 1), index_data(iter_px, 2)) = medianVecDouble_hpp(Mvoisin);}
        
         }  
   
   return(wrap(Mres));
}

