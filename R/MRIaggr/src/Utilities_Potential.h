// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> 
#include <iostream> 
#include <Rmath.h> 
#include "Utilities_General.h"
#include "Utilities_W.h"

using namespace Rcpp ;
using namespace std ;
using namespace arma ;

// fct 1 : inline List calcMultiPotential_hpp(const S4& W_SR, const S4& W_LR, const vector < double>& sample, 
// double threshold, const arma::mat& coords, 
// vector < double > distance_ref, int nbGroup_min, bool multiV, 
// double neutre);
inline vector < double > calcPotential_hpp(const vector < int>& W_LRi, const vector < int>& W_LRp, const vector < double>& W_LRx, 
                                        const vector < double>& sample, const vector < int>& subset, int n_subset, double rayon_Potentiel, 
                                        vector < double > distance_ref, double neutre);

#ifndef __VUTILITIES__
#define __VUTILITIES__

//  fct 1 ////////////////////////////////////////////////////////////
inline List calcMultiPotential_hpp(const S4& W_SR, const S4& W_LR, const vector < double>& sample, 
                                   double threshold, const arma::mat& coords, 
                                   vector < double > distance_ref, int nbGroup_min, bool multiV, 
                                   double neutre){ 
  
  int n = sample.size();
  vector < int > W_LRi = W_LR.slot("i");
  vector < int > W_LRp = W_LR.slot("p");
  vector < double > W_LRx = W_LR.slot("x");
  double maxDist = 0;
  for(int iter_x = 0 ; iter_x < n ; iter_x++){ // find maximum distance but should be pass in argument
    if(maxDist < W_LRx[iter_x]){maxDist = W_LRx[iter_x];}
  }
  
  vector < double > V(n);
  std::fill(V.begin(), V.end(), 0);
  
  //// identification des groupes spatiaux
  int n_groups ; 
  vector < int > Groupes(n);
  
  if(multiV){   
    vector < int > size_groups(0) ; 
    vector < int > size_groups_with0;
    
    vector < int > Vsubset(0) ; // ne conserve que les voxels au dessus du seuil 
    for(int iter_px = 0 ; iter_px < n ; iter_px++){
      if(sample[iter_px] >= threshold){
        Vsubset.push_back(iter_px);
      }    
    }
    
    List resR = calcGroupsW_hpp( Rcpp::as < std::vector < int > >(W_SR.slot("i")), 
                                 Rcpp::as < std::vector < int > >(W_SR.slot("p")), 
                                 Vsubset, 10000, false);
    
    Groupes = resR[0];   
    size_groups_with0 = resR[2];
    n_groups = resR[3];   
    
    //// elimination des petits groupes  spatiaux
    vector < int > Vconvert(n_groups, -1);
    int count_group = 0;
    
    for(int iter_group = 0 ; iter_group < n_groups ; iter_group++){ // maj taille groupe
      if(size_groups_with0[iter_group] >= nbGroup_min){
        count_group++;
        size_groups.push_back( size_groups_with0[iter_group] );  
        Vconvert[iter_group] = count_group;
      }else{
        size_groups_with0[iter_group] = 0;
      }    
    }
    n_groups = size_groups.size();
    
    for(int iter_px = 0 ; iter_px < n ; iter_px++){   // maj group
      if(Groupes[iter_px] > 0){
        if(size_groups_with0[Groupes[iter_px] - 1] > 0){
          Groupes[iter_px] = Vconvert[Groupes[iter_px] - 1];
        }else{
          Groupes[iter_px] = -1;
          V[iter_px] = neutre;
        }
      }
    }

  }else{
    vector < int > size_groups(0) ; 
    size_groups.push_back(0);
    
    for(int iter_px = 0 ; iter_px < n ; iter_px++){     
      if(sample[iter_px] >= threshold){
        Groupes[iter_px] = 1 ;
        size_groups[0] += 1 ;
      }else{
        Groupes[iter_px] = -1 ;
      }
    }  
    n_groups = 1;  
  }
  
  //// evaluation du rayon de chaque groupe
  vector < double > group_radius(n_groups);
  vector < bool > subset_group(n);
    
    for(int iter_g = 0 ; iter_g < n_groups ; iter_g++){
      for(size_t iter_px = 0 ; iter_px < subset_group.size(); iter_px++){
        subset_group[iter_px] = Groupes[iter_px] == (iter_g + 1);
      }
      
      group_radius[iter_g] = calcRadiusFast_hpp(coords, sample, threshold, subset_group);
    }
  
  
  //// association des pixels au groupe le plus proche
  arma::mat MppGroupes(n, n_groups); // indicatrice de la pertinence de chaque voxel dans le calcul du potentiel de chaque groupe spatial
  std::fill(MppGroupes.begin(), MppGroupes.end(), 0);  
  
  vector < double > size_groupsV(n_groups); // taille de chaque sous groupe de potentiel
  std::fill(size_groupsV.begin(), size_groupsV.end(), 0);
  
  double dist_vois;
    int groupe_vois;
    
    for(int iter_px = 0; iter_px < n ; iter_px++){
      
      for(int iter_vois = W_LRp[iter_px] ; iter_vois < W_LRp[iter_px + 1] ; iter_vois++){
       
        groupe_vois = Groupes[W_LRi[iter_vois]] - 1; // groupe spatial de lesion  ou W_LRi[iter_vois] indique l indice de iter_vois, le -1 vient du fait que le groupes sont numerotes a partir de 1
        
        if(groupe_vois>= 0 && MppGroupes(iter_px, groupe_vois) == 0){ // -2 indique un voxel en dessous du seuil dont la contribution au potentiel est negligeable
          
          dist_vois = distance_ref[W_LRx[iter_vois]]; 
         
          if(dist_vois <= group_radius[groupe_vois]){ 
            MppGroupes(iter_px, groupe_vois) = 1 ;        
            size_groupsV[groupe_vois]++ ;       
          }
        }
      }
      
    }
  
  
  // evaluation du potentiel pour les differents groupes
  {vector < double > V_tempo;
   vector < int > subset_g;
    int position, index_px;
    
    for(int iter_g = 0; iter_g < n_groups; iter_g++){
      
      subset_g.resize(size_groupsV[iter_g]);
      position = 0;
      
      // vecteur indiquant le groupe a considerer
      for(int iter_px = 0; iter_px < n ; iter_px++){
        if(MppGroupes(iter_px, iter_g) == 1){
          subset_g[position] = iter_px; 
          position++;
        }
      }
      
      V_tempo = calcPotential_hpp(W_LRi, W_LRp, W_LRx, sample, subset_g, size_groupsV[iter_g], group_radius[iter_g], distance_ref, neutre);
      
      // attribution du potentiel le plus eleve
      for(size_t iter_px = 0; iter_px < subset_g.size() ; iter_px++){
        index_px = subset_g[iter_px];
        if(V_tempo[iter_px] > V[index_px])
        {V[index_px] = V_tempo[iter_px];}
      }      
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("V") = V, 
                            Rcpp::Named("Radius") = group_radius);
}

//  fct 2 ////////////////////////////////////////////////////////////
inline vector < double > calcPotential_hpp(const vector < int>& W_LRi, const vector < int>& W_LRp, const vector < double>& W_LRx, 
                                        const vector < double>& sample, const vector < int>& subset, int n_subset, double rayon_Potentiel, 
                                        vector < double > distance_ref, double neutre){
  
  // preparation  
  vector < double > V(n_subset); 
  
  // cas particulier 
  if(n_subset == 1){
    std::fill(V.begin(), V.end(), neutre);
    return V;
  }
  
  vector < double > distance;
  vector < double > distance_min;
  
  // distances retenues
  unsigned int iter = 0;
  bool criteria = true;
  
  // ordre de voisinage a considerer
  while(criteria && iter < distance_ref.size()){       
    criteria = distance_ref[iter] <= rayon_Potentiel;
    if(criteria)
    {distance.push_back(distance_ref[iter]);
      if(iter == 0)
      {distance_min.push_back(-0.1);}else
      {distance_min.push_back(distance_ref[iter - 1]);}
    }
    iter++;
  }
  int n_distance = distance.size();
  
  // cas pathologique
  if(n_distance == 0){
    std::fill(V.begin(), V.end(), neutre);
    return V ;
  }
  
  // evaluation de la densite multi echelle 
  vector < double > densite(n_distance);  
  vector < int > nb_Vois(n_distance);
  int index_px;
  double norm;
  
    for(int iter_px = 0; iter_px < n_subset; iter_px++){
      norm = 0;
      std::fill(densite.begin(), densite.end(), 0); 
      std::fill(nb_Vois.begin(), nb_Vois.end(), 0);
      index_px = subset[iter_px];
      
      for(int iter_vois = W_LRp[index_px]; iter_vois < W_LRp[index_px + 1]; iter_vois++){
        
        if(W_LRx[iter_vois] < n_distance){// && index_px != W_LRi[iter_vois]){   // contribution de chaque voisin aux densites
          nb_Vois[W_LRx[iter_vois]]++;
          densite[W_LRx[iter_vois]]+= sample[W_LRi[iter_vois]]; 
        }
      }
      
      for(int iter_d = 0; iter_d < n_distance; iter_d++){   // evaluation des densites
        if(nb_Vois[iter_d] > 0){        
          densite[iter_d]/=nb_Vois[iter_d]; //
          V[iter_px] += densite[iter_d];
          norm++;        
        }
      }
      
      V[iter_px]/=norm;
    }
  
  
  return V;
}

#endif //__VUTILITIES__
