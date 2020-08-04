// [[Rcpp::depends(RcppArmadillo,RcppProgress)]]
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <iostream>
#include <Rmath.h>
#include "Utilities_General.hpp" // pour CST_PROGRESS

using namespace Rcpp ;
using namespace std ;
using namespace arma ;



// fct 1 : inline List calcGroupsCoords_hpp(const arma::mat& coords_NNA, const vector<int>& index_NNA, const arma::mat& Neighborhood, vector<int> coords_max,
//                           int max_groups, bool trace);
// fct 2 : inline List calcGroupsW_hpp(const S4& W, const vector<int>& subset, int max_groups);
inline  List subsetW_hpp(const vector<int>& W_i, const vector<int>& W_p, vector<int> subset); // fct 3 :
// fct 4 : inline List calcRadius_hpp(const arma::mat& coords, const vector<double>& sample, double threshold, const vector<bool>& subset_bary, bool trace);
// fct 4bis : inline double calcRadiusFast_hpp(const arma::mat& coords, const vector<double>& sample, double threshold, const vector<bool>& subset_bary);
// List calcBlockW_hpp(const vector<int>& W_i, const vector<int>& W_p, const vector<int>& site_order, bool trace)

#ifndef __WUTILITIES__
#define __WUTILITIES__

//  fct 1 ////////////////////////////////////////////////////////////
inline List calcGroupsCoords_hpp(const arma::mat& coords_NNA, const vector<int>& index_NNA, const arma::mat& Neighborhood, vector<int> coords_max,
                                 int max_groups, bool trace){
  
  // trace function
  Progress testUser(trace*CST_PROGRESS, trace);
  double value_trace = 0 ;
  
  //
  int n=coords_NNA.n_rows; // nb de vx
  int p_Mneighbors=Neighborhood.n_cols; // dim du voisinage
  int n_Mneighbors=Neighborhood.n_rows; // taille du voisinage
  vector< int> cum_coords_max(coords_max.size());
  cum_coords_max[0] = 1;
  for (size_t iter=1; iter<coords_max.size(); iter++){cum_coords_max[iter]=cum_coords_max[iter-1]*coords_max[iter-1];}  
  
  int iter_group=0; // iteration sur les groupes
  int iter_res; // iteration pour le stockage
  int n_neighbors,n_neighbors_new; // nb de nouveaux voisins
  int Sum_group_size=0 ;
  int index_tempo ;
  double sum_tempo ; // int sum_tempo ;
  size_t iter_vois_min ;
  
  vector<int> group_size(0); // taille de chaque groupe
  vector<int> residual(n); // index des voxels a classer
  vector<int> coords_newVoisin(0);
  for(int iter=0 ; iter<n ; iter++){residual[iter] = iter;}
  int n_allocated=0; // iteration sur les voxels
  arma::mat group(n,2);  // definition des groupes  
  bool test_newVoisin;
  
  while(n_allocated<n && iter_group<=max_groups){ // tant qu il reste des vx a classer et que le nb max de groupe n est pas atteint
    iter_group = iter_group + 1 ;
    
    // ajout d un nouveau groupe
    group(n_allocated,0) = residual[0];
    group(n_allocated,1) = iter_group ;
    n_allocated ++;  
    residual.erase(residual.begin()); // elimine le voxel affecte des voxels restants    
    
    group_size.push_back(group_size.size());
    group_size[group_size.size()-1] = 1 ;
    n_neighbors = 1;
    
    while(n_neighbors>0 && n_allocated<n){ // tant qu il reste des nouveaux voisins
      n_neighbors_new=0;
      coords_newVoisin.resize(0);
      
      for(int iter_vois=Sum_group_size+group_size[iter_group-1]-n_neighbors ; iter_vois < Sum_group_size+group_size[iter_group-1] ; iter_vois++){ // pour chaque voisin
        
        for(int iter_l=0 ; iter_l < n_Mneighbors ; iter_l++){ // identifier les nouveaux voisins potentiels
          test_newVoisin = true;    
          index_tempo=0;
          
          for(int iter_p=0 ; iter_p < p_Mneighbors ; iter_p++){
            
            sum_tempo = coords_NNA(group(iter_vois,0),iter_p) + Neighborhood(iter_l,iter_p);
            
            if(sum_tempo <0 || sum_tempo >= coords_max[iter_p] ){ // rejet immediat si hors du domaine
              test_newVoisin = false;
              break;                     
            }else{
              index_tempo += sum_tempo*cum_coords_max[iter_p] ;
            }                 
            
          }
          
          if(test_newVoisin){   // indice du nouveau voising potentiel
            coords_newVoisin.push_back(coords_newVoisin.size());  
            coords_newVoisin[coords_newVoisin.size()-1] = index_tempo;
          } 
        }
      }  
      
      std::sort(coords_newVoisin.begin(), coords_newVoisin.end());
      iter_vois_min=0;       
      
      iter_res = 0 ;
      while(iter_res<n-n_allocated && iter_vois_min<coords_newVoisin.size()){
        
        while(index_NNA[residual[iter_res]]>coords_newVoisin[iter_vois_min] && iter_vois_min<coords_newVoisin.size()){
          coords_newVoisin.erase(coords_newVoisin.begin() + iter_vois_min); // iter_vois_min++;                
        }
        
        if(index_NNA[residual[iter_res]]==coords_newVoisin[iter_vois_min]){ // matching               
          // attribution du groupe au voxel
          group(n_allocated,0) = residual[iter_res];
          group(n_allocated,1) = iter_group ;
          n_allocated++ ;  
          
          // maj de residual et neighbors_new
          residual.erase(residual.begin() + iter_res);
          n_neighbors_new++ ;
          iter_vois_min++;                                
        }else{
          iter_res++;
        }              
      }
      
      n_neighbors = n_neighbors_new;
      group_size[iter_group-1] += n_neighbors;
      
      if(trace && n_allocated>=value_trace){
        value_trace = min(1.0*n,value_trace + n/CST_PROGRESS);
        testUser.increment();
      }
      if (Progress::check_abort() ){
        return Rcpp::List::create(Rcpp::Named("group") = NA_REAL,
                                  Rcpp::Named("group_size") = NA_REAL
        );
      }
      
    }
    
    Sum_group_size += group_size[iter_group-1];
    
  }  

  //  export
  return(List::create(Named("group")  = group,
                      Named("group_size")  = group_size)                     
  );
  
}

// fct 2 ////////////////////////////////////////////////////////////
inline List calcGroupsW_hpp(const vector<int>& W0_i, const vector<int>& W0_p, const vector<int>& subset, int max_groups){
  
  int n;
  
  if(subset[0]>=0){
    n = subset.size(); 
  }else{
    n = W0_p.size()-1;
  }
  
  // cas pathologique
  if(n==1){    
    return(List::create(Named("group")  = 1,
                        Named("group_subset")  = NA_REAL,
                        Named("group_size")  = 1,
                        Named("nb_groups") = 1,
                        Named("n")  = 1)                     
    );
  }
  
  // subset
  int n_all;
  vector<int> W_i;
  vector<int> W_p;
  
  if(subset[0]>=0){
    n_all = W0_p.size()-1;
    
    List res_subset=subsetW_hpp(W0_i,W0_p,subset);
    
    W_i = res_subset[0];  
    W_p = res_subset[1]; 
  }else{
    W_i = W0_i;  // \Achanger car pas optimal du tout
    W_p = W0_p;  // \Achanger car pas optimal du tout
  }
  
  // initialization
  vector<int> group(n,-1); 
  vector<int> groupeSize(0);  
  vector<int> indexObs(n); 
  
  // loop
  for(int iter=0 ; iter<n ; iter++){indexObs[iter]=iter;}
  
  vector<int> voisin_0; 
  vector<int> voisin_1; 
  
  size_t iter_newvois, iter_index;
  
  // iteration sur les groupes
  int iter_group=0; 
  
  while(iter_group < max_groups && indexObs.size()>0){
    
    group[indexObs[0]] = iter_group + 1 ; 
    groupeSize.push_back(1);
    
    voisin_0.resize(1);
    voisin_0[0]=indexObs[0]; 
    indexObs.erase( indexObs.begin() );
    
    while(voisin_0.size()>0 && indexObs.size()>0){ // tant qu il y a des nouveaux voisins
      voisin_1.resize(0);
      
      for(size_t iter_vois=0 ; iter_vois<voisin_0.size() ; iter_vois++){ // chercher les nouveaux voisins
        
        if(abs(W_p[voisin_0[iter_vois]] - W_p[voisin_0[iter_vois]+1])>0.1){ // si le point a des voisins
          
          for(int iter_newvois=W_p[voisin_0[iter_vois]]; iter_newvois<W_p[voisin_0[iter_vois]+1]; iter_newvois++){
            voisin_1.push_back(W_i[iter_newvois]);  
          }
        }
        
      }    
      
      std::sort(voisin_1.begin(), voisin_1.end()); // rangement par ordre croissant
      
      iter_newvois=0;
      iter_index=0;
      
      while(iter_newvois<voisin_1.size() && iter_index < indexObs.size()){
        
        while(voisin_1[iter_newvois] < indexObs[iter_index] && iter_newvois<voisin_1.size()){
          voisin_1.erase(voisin_1.begin() + iter_newvois); 
        }     
        
        if(voisin_1[iter_newvois] == indexObs[iter_index]){ 
          
          group[voisin_1[iter_newvois]] = iter_group + 1;
          indexObs.erase(indexObs.begin() + iter_index);   
          iter_newvois++; 
          
        }else{
          iter_index++;
        }
        
      }
      
      if(iter_newvois<voisin_1.size()){    
        voisin_1.erase(voisin_1.begin() + iter_newvois,voisin_1.begin() + voisin_1.size());          
      }
      
      voisin_0 = voisin_1;
      groupeSize[iter_group] += voisin_0.size();
    }
    
    iter_group++;
    
  }
  
  // export
  int group_length=groupeSize.size();
  int nbObsRes=indexObs.size();
  
  if(subset[0]<0){
    return(List::create(
        Named("group")  = group,
        Named("group_subset")  = NA_REAL,                    
        Named("group_size")  = groupeSize,
        Named("nb_groups") = group_length,
        Named("n")  = nbObsRes)                     
    );
  }else{
  vector<int> group_all(n_all,-1);
    for(size_t iter_all = 0 ; iter_all < group.size() ; iter_all++){
      group_all[subset[iter_all]] = group[iter_all];
    }
    
    return(List::create(
        Named("group")  = group_all,
        Named("group_subset")  = group,                    
        Named("group_size")  = groupeSize,
        Named("nb_groups") = group_length,
        Named("n")  = nbObsRes)                     
    );
  }
   
}

// fct 3 ////////////////////////////////////////////////////////////
inline List subsetW_hpp(const vector<int>& W_i, const vector<int>& W_p, vector<int> subset){
  
  std::sort(subset.begin(), subset.end());
  
  int n_subset = subset.size();
  int index_min,index_max;
  vector<int> correspondance(W_p.size()-1,0);
  for(int iter_subset=0; iter_subset<n_subset ; iter_subset++){
    correspondance[subset[iter_subset]] = iter_subset;
  }
  vector<int> index_i,index_i_sort;
  vector<int> Wnew_p(n_subset+1,0);
  vector<int> Wnew_i;  
  
  for(int iter_subset=0; iter_subset<n_subset ; iter_subset++){
    
    // ne regarde que les observations d interets
    index_min = W_p[subset[iter_subset]];
    index_max = W_p[subset[iter_subset]+1];
    
    index_i.resize(index_max-index_min);
    for(int iter_i=index_min ; iter_i<index_max ; iter_i++){
      index_i[iter_i-index_min]= W_i[iter_i];
    }
    
    index_i_sort.resize(0);
    std::sort(index_i.begin(), index_i.end());
    
    // restreint les voisins aux observations d interet
    set_intersection(index_i.begin(),index_i.end(),subset.begin(),subset.end(),back_inserter(index_i_sort));
    
    if(index_i_sort.size()>0){
      for(size_t iter_i_sort=0 ; iter_i_sort < index_i_sort.size() ; iter_i_sort++){
        Wnew_i.push_back(correspondance[index_i_sort[iter_i_sort]]);
        index_i_sort[iter_i_sort] = correspondance[index_i_sort[iter_i_sort]];
      }
    }
    
    Wnew_p[iter_subset+1] = Wnew_p[iter_subset]+index_i_sort.size();
    
  }
  
  return(List::create(Named("W_i")  = Wnew_i,
                      Named("W_p")  = Wnew_p));       
  
}

// fct 4 ////////////////////////////////////////////////////////////
inline List calcRadius_hpp(const arma::mat& coords, const vector<double>& sample, double threshold, const vector<bool>& subset_bary, bool trace){
  
  // preparation
  int n = sample.size();
  vector<int> Groupe_potentiel; 
  for(int iter_n=0; iter_n<n; iter_n++)
  { if(sample[iter_n]>threshold && subset_bary[iter_n])
  {Groupe_potentiel.push_back(iter_n);}       
  }
  int n_Gp = Groupe_potentiel.size();
  
  // cas particulier
  if(n_Gp<2){
    if(trace==true){Rcout << "Rayon : " << 0 << endl;}  
    return(0);
  }
  
  int index_Gp;     
  double Rayon_medAVC=0, norm=0;
  
  // calcul du rayon - moyenne ponderee des distances au barycentre
  int D = coords.n_cols;
  
  vector<double> coord_bary(D);
  std::fill(coord_bary.begin(),coord_bary.end(),0);  
  // calcul du barycentre
  for(int iter_Gp=0 ; iter_Gp<n_Gp ; iter_Gp++){ 
    index_Gp = Groupe_potentiel[iter_Gp];  
    
    for( int iter_d=0; iter_d < D ; iter_d++){
      coord_bary[iter_d] += coords(index_Gp,iter_d)*sample[index_Gp];
    }
    norm += sample[index_Gp];
  }
  for( int iter_d=0; iter_d < D ; iter_d++){
    coord_bary[iter_d] = coord_bary[iter_d]/norm;
  }
  
  // moyenne des distances
  norm=0; 
  vector<double> distance(n,0);
  for(int iter_Gp=0 ; iter_Gp<n_Gp ; iter_Gp++){ 
    index_Gp = Groupe_potentiel[iter_Gp];
    
    for( int iter_d=0; iter_d < D ; iter_d++){
      distance[index_Gp] += pow(coord_bary[iter_d]-coords(index_Gp,iter_d),2);      
    }
    distance[index_Gp] = sqrt(distance[index_Gp]);
    Rayon_medAVC += distance[index_Gp]*sample[index_Gp];
    norm += sample[index_Gp];
  }
  
  
  Rayon_medAVC = Rayon_medAVC/norm ;
  if(trace==true){Rcout << "Rayon : " << Rayon_medAVC << endl;}
  
  return(List::create(Named("radius") = Rayon_medAVC,
                      Named("distance") = distance,
                      Named("coord_bary") = coord_bary));   
}

// fct 4bis ////////////////////////////////////////////////////////////
inline double calcRadiusFast_hpp(const arma::mat& coords, const vector<double>& sample, double threshold, const vector<bool>& subset_bary){
  
  // preparation
  int n = sample.size();
  vector<int> Groupe_potentiel; 
  for(int iter_n=0; iter_n<n; iter_n++){
    if(sample[iter_n]>threshold && subset_bary[iter_n])
    {Groupe_potentiel.push_back(iter_n);}       
  }
  int n_Gp = Groupe_potentiel.size();
  
  // cas particulier
  if(n_Gp<2){
    return(0);
  }
  
  int index_Gp;     
  double Rayon_medAVC=0, norm=0;
  
  // calcul du rayon - moyenne ponderee des distances au barycentre
  int D = coords.n_cols;
  
  vector<double> coord_bary(D);
  std::fill(coord_bary.begin(),coord_bary.end(),0);  
  // calcul du barycentre
  for(int iter_Gp=0 ; iter_Gp<n_Gp ; iter_Gp++){ 
    index_Gp = Groupe_potentiel[iter_Gp];  
    
    for( int iter_d=0; iter_d < D ; iter_d++){
      coord_bary[iter_d] += coords(index_Gp,iter_d)*sample[index_Gp];
    }
    norm += sample[index_Gp];
  }
  for( int iter_d=0; iter_d < D ; iter_d++){
    coord_bary[iter_d] = coord_bary[iter_d]/norm;
  }
  
  // moyenne des distances
  norm=0;
  
  double distance;
  for(int iter_Gp=0 ; iter_Gp<n_Gp ; iter_Gp++){ 
    index_Gp = Groupe_potentiel[iter_Gp];
    distance = 0 ;
    
    for( int iter_d=0; iter_d < D ; iter_d++){
      distance += pow(coord_bary[iter_d]-coords(index_Gp,iter_d),2);      
    }
    Rayon_medAVC += sqrt(distance)*sample[index_Gp];
    norm += sample[index_Gp];
  }
  Rayon_medAVC = Rayon_medAVC/norm ;
  
  
  return(Rayon_medAVC);
}

// fct 5 ////////////////////////////////////////////////////////////
inline List calcBlockW_hpp(const vector<int>& W_i, const vector<int>& W_p, const vector<int>& site_order, bool trace){
  
  // trace function
  Progress testUser(trace*CST_PROGRESS, trace);
  double value_trace = 0 ;
  
  //
  int n = W_p.size()-1; 
  int n_groups;
  int n_neighbors;
  int site;
  int iter_group;
  vector<vector<int> > ls_groups;
  vector<int> size_groups;
  vector<int> sites;
  bool test_neighbor ;
  bool no_site_order = (site_order[0]<0);
  
  // initialisation
  ls_groups.push_back(sites);
  ls_groups[0].push_back(0);
  size_groups.push_back(1);
  n_groups=1 ;
  
  // main loop
  for(int iter_site=1 ; iter_site<n ; iter_site++){
    
    if(no_site_order){site = iter_site;}else{site = site_order[iter_site];}
    
    n_neighbors = W_p[site+1]-W_p[site] ;
    
    test_neighbor = true ; // just to launch the loop
    iter_group = -1;
    
    while(test_neighbor && iter_group<(n_groups-1)){
      iter_group++;
      test_neighbor = false;
      
      for(int iter_neighbors=0 ; iter_neighbors<n_neighbors; iter_neighbors++){
        for(int iter_sitegroup=0 ; iter_sitegroup<size_groups[iter_group] ; iter_sitegroup++){
          if(abs(ls_groups[iter_group][iter_sitegroup]-W_i[W_p[site]+iter_neighbors])<0.5){test_neighbor = true;}
        }
      } 
    }
    
    if(test_neighbor==false){ // add site to existing group
      ls_groups[iter_group].push_back(site);
      size_groups[iter_group]++;
    }else{ // create a new group
      ls_groups.push_back(sites);
      ls_groups[n_groups].push_back(site);
      size_groups.push_back(1);
      n_groups++;
    }
    
    if(trace && iter_site>=value_trace){
      value_trace = min(1.0*n,value_trace + n/CST_PROGRESS);
      testUser.increment();
    }
    if (Progress::check_abort() ){
      return Rcpp::List::create(Rcpp::Named("ls_groups") = NA_REAL,
                                Rcpp::Named("size_groups") = NA_REAL,
                                Rcpp::Named("n_groups") = NA_REAL
      );
    }
    
  }
  
  return Rcpp::List::create(Rcpp::Named("ls_groups") = ls_groups,
                            Rcpp::Named("size_groups") = size_groups,
                            Rcpp::Named("n_groups") = n_groups);    
}

#endif //__WUTILITIES__
