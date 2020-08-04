// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> 
#include <iostream> 
#include <Rmath.h> 
#include "Utilities_General.h"


using namespace Rcpp ;
using namespace std ;
using namespace arma ;

// List calcHemi_cpp(const IntegerVector& coordsI, const IntegerVector& coordsJ, List ls_indexK, int n_num, const NumericVector& value, int n, 
//                  double i_pos, double j_pos, double angle_pos, 
//                  double penaltyNA, double sd_data, double p, bool symetrie)
// List calcContro_cpp(const arma::mat& contrast, const arma::mat& coords_px, const IntegerVector& index_k, const IntegerVector& index_k_contro, 
//                    double d_lim, double lambda, int param_ref, double var_ref , 
//                    bool type_moy, bool type_med, bool type_NN, bool verbose)

//  1 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List calcHemi_cpp(const IntegerVector& coordsI, const IntegerVector& coordsJ, List ls_indexK, int n_num, const NumericVector& value, int n, 
                  double i_pos, double j_pos, double angle_pos, 
                  double penaltyNA, double sd_data, double p, bool symetrie){ 
  
  // validite de angle_pos  
  if(abs(angle_pos) > M_PI / 2){
    if(symetrie){
      return(List::create(Named("criteria")  = NA_REAL, 
                          Named("numberAssociated")  = NA_REAL, 
                          Named("pcNA")  = NA_REAL,                      
                          Named("compromise")  = -1));
      
    }else{
      return(List::create(Named("criteria")  = NA_REAL, 
                          Named("numberAssociated")  = NA_REAL, 
                          Named("pcNA")  = NA_REAL,                      
                          Named("compromise")  = + 100000));
    }
  }
  
  // critere de symetrie
  double valeur_vois, dist_vois, asymetrie_valeur = 0, poids_vois;
  int asymetrie_nb = 0, voxel_nb = 0;
  double distI, distJ ;
  double sdp_data = pow(sd_data, p);
  
  // parametrage du nouveau de repere
  int n_indexK, n_hemiL, n_hemiR, iter_vxAll ;
  double cos_alpha = cos(angle_pos);
  double sin_alpha = sin(angle_pos);
  
  vector < int > indexK;
  vector < double > i_nouveauL, i_nouveauR, j_nouveauL, j_nouveauR, valueL, valueR;
  double i_tempo;
  
  for(int iter_k = 0 ; iter_k < n_num ; iter_k++){
    
    indexK = ls_indexK[iter_k];
    n_indexK = indexK.size();
    n_hemiL = 0;
    n_hemiR = 0;
    
    i_nouveauL.resize(n_indexK);
    i_nouveauR.resize(n_indexK);
    j_nouveauL.resize(n_indexK);
    j_nouveauR.resize(n_indexK);
    valueL.resize(n_indexK);
    valueR.resize(n_indexK);
    
    for(int iter_vx = 0 ; iter_vx < n_indexK ; iter_vx++){
      iter_vxAll = indexK[iter_vx];  
      
      // nouvelles coordonnees
      distI = coordsI[iter_vxAll]-i_pos;
      distJ = coordsJ[iter_vxAll]-j_pos;
      i_tempo = cos_alpha * distI - sin_alpha * distJ;
      
      if(i_tempo > 0){ // left hemisphere (radiological convention)
      
      i_nouveauL[n_hemiL] = i_tempo;
      j_nouveauL[n_hemiL] = sin_alpha * distI + cos_alpha * distJ;
      valueL[n_hemiL] = value[iter_vxAll];
      n_hemiL++;
      
      }else{ // right hemisphere
      
      i_nouveauR[n_hemiR] = i_tempo;
      j_nouveauR[n_hemiR] = sin_alpha * distI + cos_alpha * distJ;
      valueR[n_hemiR] = value[iter_vxAll];      
      n_hemiR++;
      
      }
      
    }
    
    // evalutation du critere de symetrie vis a vis de l hemisphere gauche
    voxel_nb += n_hemiL;
 
    for(int iter_pxL = 0 ; iter_pxL < n_hemiL ; iter_pxL++){
      poids_vois = 0;
      valeur_vois = 0;        
            
      for(int iter_pxR = 0 ; iter_pxR < n_hemiR ; iter_pxR++){ // iteration sur les voisins : moyenne ponderee par les distance des valeurs voisines
        
        distI = abs(i_nouveauL[iter_pxL]+i_nouveauR[iter_pxR]) ;
        distJ = abs(j_nouveauL[iter_pxL]-j_nouveauR[iter_pxR]) ;
        
        if(distI < 1 && distJ < 1){
          dist_vois = sqrt(pow(distI, 2) + pow(distJ, 2));
          valeur_vois += valueR[iter_pxR] * (2 - dist_vois)         ;           
          poids_vois += (2 - dist_vois) ;
        }
        
      }
      
      if(poids_vois > 0){
        
        valeur_vois = valeur_vois / poids_vois;
        if(symetrie){
          asymetrie_valeur +=  1 / (1 + pow(abs(valueL[iter_pxL]-valeur_vois), p) / sdp_data);
        }else{        
          asymetrie_valeur += pow(abs(valueL[iter_pxL]-valeur_vois), p) / sdp_data   ;
        }
         
        asymetrie_nb++;
      }
    }
    
    // evalutation du critere de symetrie vis a vis de l hemisphere droit
    voxel_nb += n_hemiR;
 
    for(int iter_pxR = 0 ; iter_pxR < n_hemiR ; iter_pxR++){
      poids_vois = 0;
      valeur_vois = 0;        
            
      for(int iter_pxL = 0 ; iter_pxL < n_hemiL ; iter_pxL++){ // iteration sur les voisins : moyenne ponderee par les distance des valeurs voisines
        
        distI = abs(i_nouveauR[iter_pxL] + i_nouveauL[iter_pxR]) ;
        distJ = abs(j_nouveauR[iter_pxL] - j_nouveauL[iter_pxR]) ;
        
        if(distI < 1 && distJ < 1){
          dist_vois = sqrt(pow(distI, 2) + pow(distJ, 2));
          valeur_vois += valueL[iter_pxL] * (2 - dist_vois)         ;           
          poids_vois += (2 - dist_vois) ;
        }
        
      }
      
      if(poids_vois > 0){
        
        valeur_vois = valeur_vois / poids_vois;
        if(symetrie){
          asymetrie_valeur += 1 / (1 + pow(abs(valueR[iter_pxR]-valeur_vois), p) / sdp_data);
        }else{        
          asymetrie_valeur += pow(abs(valueR[iter_pxR]-valeur_vois), p) / sdp_data   ;
        }
         
        asymetrie_nb++;
      }
    }
    
  }
  
  // no voxel 
  if(asymetrie_nb == 0){
    if(symetrie){
      return(List::create(Named("criteria") = NA_REAL, 
                          Named("numberAssociated") = 0, 
                          Named("pcNA") = 1,                      
                          Named("compromise") = -1));
      
    }else{
      return(List::create(Named("criteria") = NA_REAL, 
                          Named("numberAssociated") = 0, 
                          Named("pcNA") = 1,                      
                          Named("compromise") =  +100000));
    }
  }
  
  // bilan
  double compromise;  
  double pcNA = (voxel_nb-asymetrie_nb) / (double)voxel_nb;
  if(penaltyNA <= 0){
    compromise = asymetrie_valeur / asymetrie_nb;
  }else{
    if(symetrie){
      compromise = asymetrie_valeur / asymetrie_nb - pow(pcNA, penaltyNA); // NA_nb penalize for mising data
    }else{
      compromise = asymetrie_valeur / asymetrie_nb + pcNA * penaltyNA; // NA_nb penalize for mising data
    }
  }
  
  // export
  return(List::create(Named("criteria") = asymetrie_valeur, 
                      Named("numberAssociated") = asymetrie_nb, 
                      Named("pcNA") = pcNA,                      
                      Named("compromise") = compromise));
  
}


//  2 ////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List calcContro_cpp(const arma::mat& contrast, const arma::mat& coords_px, const IntegerVector& index_k, const IntegerVector& index_k_contro, 
                    double d_lim, double lambda, int param_ref, double var_ref , 
                    bool type_moy, bool type_med, bool type_NN, bool verbose){
                
                // preparation
                int n_k = index_k.size();
                int n_k_contro = index_k_contro.size();
                int D = contrast.n_cols;
                double dist1, dist2 ;
                bool voisin_proche ;
                int nb_neighbors = pow(2 * d_lim, 2);
                int index_px, index_contro, nb_vois;       
                arma::mat valeur_NormContro(n_k, D);
                std::fill(valeur_NormContro.begin(), valeur_NormContro.end(), NA_REAL);
                vector < bool > index_plot_k(n_k, false), index_plot_k_contro(n_k_contro, false); 
                vector < int > index_vois;
                vector < int > iter_vois;
                vector < double > valeurs_miroir;
                
                double poids_vois, norm , valeur_moyenne ; 
                int index_ref_contro; 
                double diff_ref, diff_tempo, dist_vois;     
                  
                // iteration sur les pixels
                for(int iter_px = 0 ; iter_px < n_k ; iter_px++){
                  index_px = index_k[iter_px];
                  index_vois.resize(0); index_vois.reserve(nb_neighbors);
                  iter_vois.resize(0); iter_vois.reserve(nb_neighbors);
                  voisin_proche = false;
                                
                  for(int iter_contro = 0 ; iter_contro < n_k_contro ; iter_contro++ ){
                  index_contro = index_k_contro[iter_contro];
                  // evaluation des distances
                  dist1 = abs(coords_px(index_contro, 0) + coords_px(index_px, 0)) ;
                  dist2 = abs(coords_px(index_contro, 1) - coords_px(index_px, 1)) ;
                  
                 
                  // si le point contro est assez proche
                  if( dist1 <= d_lim && dist2 <= d_lim )
                  {   
                      if(sqrt(pow(dist1, 2) + pow(dist2, 2)) <= 1)
                     {voisin_proche = true ;}
                  
                    // maj des vois utilises pour la normalisation
                    index_vois.push_back(index_contro)  ;             
                    if(verbose){iter_vois.push_back(iter_contro);} 
                  }
                  
                 }
                // si l on trouve des voisins dont au moins un proche
                nb_vois = index_vois.size(); 
                if(voisin_proche && nb_vois > 0){
                
                // maj des px utilises
                if(verbose){
                index_plot_k[iter_px] = true;
                for(int iter = 0 ; iter < nb_vois ; iter++){
                index_plot_k_contro[iter_vois[iter]] = true;
                }
                }
                
                if(type_moy){
                                     
                  for(int iter_d = 0 ; iter_d < D ; iter_d++){
                    norm = 0;
                    valeur_moyenne = 0;

                    for(int iter_vois = 0 ;  iter_vois < nb_vois ; iter_vois++){
                          
                          index_contro = index_vois[iter_vois];
                          poids_vois =  2 * d_lim-sqrt(pow(coords_px(index_contro, 0) + coords_px(index_px, 0), 2) + pow(coords_px(index_contro, 1) - coords_px(index_px, 1), 2));

                          valeur_moyenne += poids_vois * contrast(index_contro, iter_d);
                          norm += poids_vois;
                    }

                     valeur_NormContro(iter_px, iter_d) = contrast(index_px, iter_d) - valeur_moyenne / norm;
                  }
                }
                
                if(type_med){
                  
                  for(int iter_d = 0 ; iter_d < D ; iter_d++){
                    
                    // vecteur des donnees
                    valeurs_miroir.resize(nb_vois);
                     for(int iter_vois = 0 ;  iter_vois < nb_vois ; iter_vois++){
                          valeurs_miroir[iter_vois] = contrast(index_vois[iter_vois], iter_d);
                     }
                    
                     valeur_NormContro(iter_px, iter_d) = contrast(index_px, iter_d) - medianVecDouble_hpp(valeurs_miroir); //medianNV_hpp(valeurs_miroir);
                  }
                }
                
                if(type_NN){
                  index_contro = index_vois[0];
                  dist_vois = sqrt(pow(coords_px(index_contro, 0) + coords_px(index_px, 0), 2) + pow(coords_px(index_contro, 1) - coords_px(index_px, 1), 2));
                  diff_ref = abs(contrast(index_px, param_ref) - contrast(index_contro, param_ref)) / var_ref + lambda * dist_vois;
                  index_ref_contro = index_contro;

                  // identification du pixel le plus proche en intensite
                  for(int iter_vois = 1 ;  iter_vois < nb_vois ; iter_vois++){
                    index_contro = index_vois[iter_vois];
                    dist_vois = sqrt(pow(coords_px(index_contro,0) + coords_px(index_px,0), 2) + pow(coords_px(index_contro, 1) - coords_px(index_px, 1), 2));
                         
                    diff_tempo = abs(contrast(index_px, param_ref) - contrast(index_contro, param_ref)) / var_ref + lambda * dist_vois;
                    
                    if(diff_tempo - diff_ref < -0.0000000000001){
                    diff_ref = diff_tempo;
                    index_ref_contro = index_contro;
                    }                  
                  }
                    for(int iter_d = 0 ; iter_d < D ; iter_d++){                  
                    valeur_NormContro(iter_px, iter_d) = contrast(index_px, iter_d) - contrast(index_ref_contro, iter_d);                 
                    }
                   
                 }
              
                }
                }
                
                return(List::create(Named("valeur_NormContro")  = valeur_NormContro, 
                          Named("index_plot_k")  = index_plot_k, 
                          Named("index_plot_k_contro")  = index_plot_k_contro
                ));
}
