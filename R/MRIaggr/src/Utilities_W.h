// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]
#include <RcppArmadillo.h> 
#include <progress.hpp> 
#include <iostream> 
#include <Rmath.h> 
#include "Utilities_General.h" // pour CST_PROGRESS

using namespace Rcpp ;
using namespace std ;
using namespace arma ;



// fct 1 : inline List calcGroupsCoords_hpp(const arma::mat& coords_NNA, const vector < int>& index_NNA, const arma::mat& Neighborhood, vector < int > coords_max, 
//                           int max_groups, bool verbose);
// fct 2 : inline List calcGroupsW_hpp(const S4& W, const vector < int>& subset, int max_groups);
inline  List subsetW_hpp(const vector < int>& W_i, const vector < int>& W_p, vector < int > subset); // fct 3 :
// fct 4 : inline List calcRadius_hpp(const arma::mat& coords, const vector < double>& sample, double threshold, const vector < bool>& subset_bary, bool verbose);
// fct 4bis : inline double calcRadiusFast_hpp(const arma::mat& coords, const vector < double>& sample, double threshold, const vector < bool>& subset_bary);
// List calcBlockW_hpp(const vector < int>& W_i, const vector < int>& W_p, const vector < int>& site_order, bool verbose)

#ifndef __WUTILITIES__
#define __WUTILITIES__


//  fct 1 ////////////////////////////////////////////////////////////
inline List calcGroupsCoords_hpp(const arma::mat& coords_NNA, vector<int> index_NNA, int min_index_NNA, int max_index_NNA,
                                 const arma::mat& Neighborhood, vector < int > coords_max, int max_groups,  bool verbose){
  
  // verbose function
  Progress testUser(verbose * CST_PROGRESS, verbose);
  double value_trace = 0 ;
  int n_allocated = 0;
  
  // preparation
  int n = coords_NNA.n_rows; // nb de vx
  int n_coords = Neighborhood.n_cols; // dim du voisinage
  int n_neighbors = Neighborhood.n_rows; // taille du voisinage
  
  vector < int > cum_coords_max(n_coords);
  cum_coords_max[0] = 1;
  for (int iter_coords = 1; iter_coords < n_coords; iter_coords++){
    cum_coords_max[iter_coords] = cum_coords_max[iter_coords - 1] * coords_max[iter_coords - 1];
  }  
  
  if(min_index_NNA > 0){
    for(int iter_obs = 0; iter_obs < n ; iter_obs++){
      index_NNA[iter_obs] -= min_index_NNA;
    }
    max_index_NNA -= min_index_NNA;
  }
  
  // pathological cases
  if(n == 0){    
    return(List::create(Named("group")  = NA_REAL, 
                        Named("group_subset")  = NA_REAL, 
                        Named("group_size")  = NA_REAL, 
                        Named("nb_groups") = 0, 
                        Named("cv")  = true)                     
    );
  }
  if(n == 1){    
    return(List::create(Named("group")  = 1, 
                        Named("group_subset")  = NA_REAL, 
                        Named("group_size")  = 1, 
                        Named("nb_groups") = 1, 
                        Named("cv")  = true)                     
    );
  }
  
  // initialisation
  vector < int > group(n, -1); 
  vector < int > groupeSize(0);  
  
  vector < bool > availableObs(max_index_NNA + 1,false);           // [TO BE THINK] may not be optimal if the field is large and the group is small
  vector < int > availableObs_indexCoords(max_index_NNA + 1, -1);  // [TO BE THINK] may not be optimal if the field is large and the group is small
  
  for(int iter_obs = 0; iter_obs < n ; iter_obs++){
    availableObs[index_NNA[iter_obs]] = true;
    availableObs_indexCoords[index_NNA[iter_obs]] = iter_obs;
  }
  
  int firstObs = 0; 
  
  vector < int > extension; 
  vector < int > extension_save; 
  
  int iter_group = 0; 
  int indexCoords_vois, coord_tempo, indexArray_newvois, indexCoord_newvois;
  bool test_newvois;
  
  while(firstObs < n && iter_group < max_groups){ // tant qu il reste des vx a classer et que le nb max de groupe n est pas atteint
    
    groupeSize.push_back(1); // add a new group
    group[firstObs] = iter_group + 1 ;  // first remaining observation is associated with the new group
    extension_save.resize(1); // new group initialized with one point
    extension_save[0] = firstObs;  // which is the first remaining observation
    availableObs[index_NNA[firstObs]] = false; // which is no more available
    
    while(firstObs < n && availableObs[index_NNA[firstObs]] == false){firstObs++;} // search for the next first remaining observation
    
    if(verbose){n_allocated += extension_save.size();}
    
    while(firstObs < n && extension_save.size() > 0){ // tant qu il reste des nouveaux voisins
      
      extension.resize(0); // potential extention for the group
      
      for(size_t iter_vois = 0 ; iter_vois < extension_save.size() ; iter_vois++){ // among the points from the previous extension
        if(firstObs >= n){break;}
        indexCoords_vois = extension_save[iter_vois]; // index of a given point from the previous extension
        
        for(int iter_newvois = 0 ; iter_newvois < n_neighbors ; iter_newvois++){ // among the neighbors of a given point from the previous extension
          if(firstObs >= n){break;}
          test_newvois = true;    
          indexArray_newvois = 0;
          
          for(int iter_coord = 0 ; iter_coord < n_coords ; iter_coord++){
            
            coord_tempo = coords_NNA(indexCoords_vois, iter_coord) + Neighborhood(iter_newvois, iter_coord);
            
            if(coord_tempo < 0 || coord_tempo >= coords_max[iter_coord] ){ // First check: observation within the domain
              test_newvois = false;
              break;                     
            }else{
              indexArray_newvois += coord_tempo * cum_coords_max[iter_coord] ;
            }                 
            
          }
          
          if(min_index_NNA>0){indexArray_newvois -= min_index_NNA;}
          
          if(test_newvois && indexArray_newvois >= 0 && indexArray_newvois <= max_index_NNA && availableObs[indexArray_newvois]){   // indice du nouveau voising potentiel dans array
            indexCoord_newvois = availableObs_indexCoords[indexArray_newvois]; 
            extension.push_back(indexCoord_newvois);  
            
            group[indexCoord_newvois] = iter_group + 1;
            availableObs[indexArray_newvois] = false;
            while(firstObs < n && availableObs[index_NNA[firstObs]] == false){firstObs++;}
          } 
        } 
      }  
      
      extension_save = extension;
      groupeSize[iter_group] += extension_save.size();
      
      // check abort
      if (Progress::check_abort() ){
        return Rcpp::List::create(Named("group")  = NA_REAL, 
                                  Named("group_subset")  = NA_REAL, 
                                  Named("group_size")  = NA_REAL, 
                                  Named("nb_groups") = NA_REAL, 
                                  Named("cv")  = false
        );
      }
      
      // display
      if(verbose){
        
        n_allocated += extension_save.size();
        
        if(n_allocated>=value_trace){
          value_trace = min(1.0 * n, value_trace + n / CST_PROGRESS);
          testUser.increment();
        }
        
      }
      
    }
    
    iter_group++;
    
  }  
  
  //  export
  return(List::create(Named("group")  = group, 
                      Named("group_subset")  = NA_REAL, 
                      Named("group_size")  = groupeSize, 
                      Named("nb_groups") = iter_group,
                      Named("cv") = firstObs >= n,
                      Named("availableObs") = availableObs_indexCoords)                     
  );
  
}

// fct 2 ////////////////////////////////////////////////////////////
inline List calcGroupsW_hpp(const vector < int>& W0_i, const vector < int>& W0_p, const vector <int>& subset, int max_groups, bool verbose){
  
  // verbose function
  Progress testUser(verbose * CST_PROGRESS, verbose);
  double value_trace = 0 ;
  int n_allocated = 0;
  
  // preparation
  int n_all, n;
  vector < int > W_i;
  vector < int > W_p;
  
  n_all = W0_p.size() - 1;
  if(subset[0] >= 0){
    n = subset.size(); 	
  }else{
    n = n_all;
  }
  
  // pathological case
  if(n == 0){    
    return(List::create(Named("group")  = NA_REAL, 
                        Named("group_subset")  = NA_REAL, 
                        Named("group_size")  = NA_REAL, 
                        Named("nb_groups") = 0, 
                        Named("cv")  = true)                     
    );
  }
  if(n == 1){    
    return(List::create(Named("group")  = 1, 
                        Named("group_subset")  = NA_REAL, 
                        Named("group_size")  = 1, 
                        Named("nb_groups") = 1, 
                        Named("cv")  = true)                     
    );
  }
  
  if(subset[0] >= 0){
    List res_subset = subsetW_hpp(W0_i, W0_p, subset);
    
    W_i = res_subset[0];  
    W_p = res_subset[1]; 
  }else{
    W_i = W0_i;  // TO BE CHANGED 
    W_p = W0_p;  // TO BE CHANGED 
  }
  
  // initialization
  vector < int > group(n, -1); 
  vector < int > groupeSize(0);  
  
  vector < bool > availableObs(n,true); 
  int firstObs = 0; 
  
  vector < int > extension; 
  vector < int > extension_save; 
  
  int newvois_tempo;
  
  // iteration sur les groupes
  int iter_group = 0; 
  
  while(firstObs < n && iter_group < max_groups){
    
    groupeSize.push_back(1); // add a new group
    group[firstObs] = iter_group + 1 ;  // first remaining observation is associated with the new group
    extension_save.resize(1); // new group initialized with one point
    extension_save[0] = firstObs;  // which is the first remaining observation
    availableObs[firstObs] = false; // which is no more available
    
    while(firstObs < n && availableObs[firstObs] == false){firstObs++;} // search for the next first remaining observation
    
    while(firstObs < n && extension_save.size() > 0){ // while the group extends and there are still remaining points, continue to extend the group
      extension.resize(0); // potential extention for the group
      
      for(size_t iter_vois = 0 ; iter_vois < extension_save.size() ; iter_vois++){ // among the points from the previous extension
        if(firstObs >= n){break;}
        
        if(abs(W_p[extension_save[iter_vois]] - W_p[extension_save[iter_vois] + 1]) > 0.1){ // if a given point from the previous extension has neighbors
          
          for(int iter_newvois = W_p[extension_save[iter_vois]]; iter_newvois < W_p[extension_save[iter_vois] + 1]; iter_newvois++){ // all of them are candidates for the new extension
            if(firstObs >= n){break;}
            
            if(availableObs[W_i[iter_newvois]]){ // if they are still available
              
              newvois_tempo = W_i[iter_newvois]; 
              
              extension.push_back(newvois_tempo);     // add them as a new extension
              group[newvois_tempo] = iter_group + 1;  // associate a group number
              
              // update remaining points
              availableObs[newvois_tempo] = false;
              
              // update first observation
              while(firstObs < n && availableObs[firstObs] == false){firstObs++;}
            }
            
          }
          
          
        }
      }  
      
      extension_save = extension;
      groupeSize[iter_group] += extension_save.size(); // update group size
      
      if (Progress::check_abort() ){
        return Rcpp::List::create(Named("group")  = NA_REAL, 
                                  Named("group_subset")  = NA_REAL, 
                                  Named("group_size")  = NA_REAL, 
                                  Named("nb_groups") = NA_REAL, 
                                  Named("cv")  = false
        );
        
        if(verbose){
          n_allocated += extension_save.size();
          
          if(n_allocated>=value_trace){
            value_trace = min(1.0 * n, value_trace + n / CST_PROGRESS);
            testUser.increment();
          }
          
        }
      }
    }
    
    iter_group++;
    
  }
  
  // export
  int group_length = groupeSize.size();
  
  if(subset[0] < 0){
    return(List::create(
        Named("group")  = group, 
        Named("group_subset")  = NA_REAL,                    
        Named("group_size")  = groupeSize, 
        Named("nb_groups") = group_length,
        Named("cv") = firstObs >= n)                     
    );
  }else{
    vector < int > group_all(n_all, -1);
    for(size_t iter_all = 0 ; iter_all < group.size() ; iter_all++){
      group_all[subset[iter_all]] = group[iter_all];
    }
    
    return(List::create(
        Named("group")  = group_all, 
        Named("group_subset")  = group,                    
        Named("group_size")  = groupeSize, 
        Named("nb_groups") = group_length,
        Named("cv") = firstObs >= n)                     
    );
  }
  
}

// fct 3 ////////////////////////////////////////////////////////////
inline List subsetW_hpp(const vector < int>& W_i, const vector < int>& W_p, vector < int > subset){
  
  std::sort(subset.begin(), subset.end());
  
  int n_subset = subset.size();
  int index_min, index_max;
  vector < int > correspondance(W_p.size() - 1, 0);
  for(int iter_subset = 0; iter_subset < n_subset ; iter_subset++){
    correspondance[subset[iter_subset]] = iter_subset;
  }
  vector < int > index_i, index_i_sort;
  vector < int > Wnew_p(n_subset + 1, 0);
  vector < int > Wnew_i;  
  
  for(int iter_subset = 0; iter_subset < n_subset ; iter_subset++){
    
    // ne regarde que les observations d interets
    index_min = W_p[subset[iter_subset]];
    index_max = W_p[subset[iter_subset] + 1];
    
    index_i.resize(index_max-index_min);
    for(int iter_i = index_min ; iter_i < index_max ; iter_i++){
      index_i[iter_i-index_min] = W_i[iter_i];
    }
    
    index_i_sort.resize(0);
    std::sort(index_i.begin(), index_i.end());
    
    // restreint les voisins aux observations d interet
    set_intersection(index_i.begin(), index_i.end(), subset.begin(), subset.end(), back_inserter(index_i_sort));
    
    if(index_i_sort.size() > 0){
      for(size_t iter_i_sort = 0 ; iter_i_sort < index_i_sort.size() ; iter_i_sort++){
        Wnew_i.push_back(correspondance[index_i_sort[iter_i_sort]]);
        index_i_sort[iter_i_sort] = correspondance[index_i_sort[iter_i_sort]];
      }
    }
    
    Wnew_p[iter_subset + 1] = Wnew_p[iter_subset]+index_i_sort.size();
    
  }
  
  return(List::create(Named("W_i")  = Wnew_i, 
                      Named("W_p")  = Wnew_p));       
  
}

// fct 4 ////////////////////////////////////////////////////////////
inline List calcRadius_hpp(const arma::mat& coords, const vector < double>& sample, double threshold, const vector < bool>& subset_bary, bool verbose){
  
  // preparation
  int n = sample.size();
  vector < int > Groupe_potentiel; 
  for(int iter_n = 0; iter_n < n; iter_n++)
  { if(sample[iter_n] > threshold && subset_bary[iter_n])
  {Groupe_potentiel.push_back(iter_n);}       
  }
  int n_Gp = Groupe_potentiel.size();
  
  // cas particulier
  if(n_Gp < 2){
    if(verbose == true){Rcout << "Rayon : " << 0 << endl;}  
    return(0);
  }
  
  int index_Gp;     
  double Rayon_medAVC = 0, norm = 0;
  
  // calcul du rayon - moyenne ponderee des distances au barycentre
  int D = coords.n_cols;
  
  vector < double > coord_bary(D);
  std::fill(coord_bary.begin(), coord_bary.end(), 0);  
  // calcul du barycentre
  for(int iter_Gp = 0 ; iter_Gp < n_Gp ; iter_Gp++){ 
    index_Gp = Groupe_potentiel[iter_Gp];  
    
    for( int iter_d = 0; iter_d < D ; iter_d++){
      coord_bary[iter_d] += coords(index_Gp, iter_d) * sample[index_Gp];
    }
    norm += sample[index_Gp];
  }
  for( int iter_d = 0; iter_d < D ; iter_d++){
    coord_bary[iter_d] = coord_bary[iter_d]/norm;
  }
  
  // moyenne des distances
  norm = 0; 
  vector < double > distance(n, 0);
  for(int iter_Gp = 0 ; iter_Gp < n_Gp ; iter_Gp++){ 
    index_Gp = Groupe_potentiel[iter_Gp];
    
    for( int iter_d = 0; iter_d < D ; iter_d++){
      distance[index_Gp] += pow(coord_bary[iter_d]-coords(index_Gp, iter_d), 2);      
    }
    distance[index_Gp] = sqrt(distance[index_Gp]);
    Rayon_medAVC += distance[index_Gp]*sample[index_Gp];
    norm += sample[index_Gp];
  }
  
  
  Rayon_medAVC = Rayon_medAVC / norm ;
  if(verbose == true){Rcout << "Rayon : " << Rayon_medAVC << endl;}
  
  return(List::create(Named("radius") = Rayon_medAVC, 
                      Named("distance") = distance, 
                      Named("coord_bary") = coord_bary));   
}

// fct 4bis ////////////////////////////////////////////////////////////
inline double calcRadiusFast_hpp(const arma::mat& coords, const vector < double>& sample, double threshold, const vector < bool>& subset_bary){
  
  // preparation
  int n = sample.size();
  vector < int > Groupe_potentiel; 
  for(int iter_n = 0; iter_n < n; iter_n++){
    if(sample[iter_n] > threshold && subset_bary[iter_n])
    {Groupe_potentiel.push_back(iter_n);}       
  }
  int n_Gp = Groupe_potentiel.size();
  
  // cas particulier
  if(n_Gp < 2){
    return(0);
  }
  
  int index_Gp;     
  double Rayon_medAVC = 0, norm = 0;
  
  // calcul du rayon - moyenne ponderee des distances au barycentre
  int D = coords.n_cols;
  
  vector < double > coord_bary(D);
  std::fill(coord_bary.begin(), coord_bary.end(), 0);  
  // calcul du barycentre
  for(int iter_Gp = 0 ; iter_Gp < n_Gp ; iter_Gp++){ 
    index_Gp = Groupe_potentiel[iter_Gp];  
    
    for( int iter_d = 0; iter_d < D ; iter_d++){
      coord_bary[iter_d] += coords(index_Gp, iter_d) * sample[index_Gp];
    }
    norm += sample[index_Gp];
  }
  for( int iter_d = 0; iter_d < D ; iter_d++){
    coord_bary[iter_d] = coord_bary[iter_d]/norm;
  }
  
  // moyenne des distances
  norm = 0;
  
  double distance;
  for(int iter_Gp = 0 ; iter_Gp < n_Gp ; iter_Gp++){ 
    index_Gp = Groupe_potentiel[iter_Gp];
    distance = 0 ;
    
    for( int iter_d = 0; iter_d < D ; iter_d++){
      distance += pow(coord_bary[iter_d]-coords(index_Gp, iter_d), 2);      
    }
    Rayon_medAVC += sqrt(distance) * sample[index_Gp];
    norm += sample[index_Gp];
  }
  Rayon_medAVC = Rayon_medAVC / norm ;
  
  
  return(Rayon_medAVC);
}

// fct 5 ////////////////////////////////////////////////////////////
inline List calcBlockW_hpp(const vector < int>& W_i, const vector < int>& W_p, const vector < int>& site_order, 
                           const vector < double>& dist_center, double dist_max, bool verbose){
  
  // verbose function
  Progress testUser(verbose * CST_PROGRESS, verbose);
  double value_trace = 0 ;
  
  //
  int n = W_p.size() - 1; 
  int n_groups;
  int n_neighbors;
  int index_neighbor_tempo;
  int site;
  int iter_group;
  double current_dist = dist_max;
  vector < vector < int > > ls_groups;
  vector < int > size_groups;
  vector < int > init_groups;
  vector < int > sites;
  bool test_neighbor ;
  bool no_site_order = (site_order[0] < 0);
  
  // initialisation
  ls_groups.push_back(sites);
  if(no_site_order){ls_groups[0].push_back(0);}else{ls_groups[0].push_back(site_order[0]);}
  size_groups.push_back(1);
  init_groups.push_back(0);
  n_groups = 1 ;
  
  // main loop
  
  for(int iter_site = 1 ; iter_site < n ; iter_site++){
    
    if(no_site_order){site = iter_site;}else{site = site_order[iter_site];}
    
    if(dist_center[site] > current_dist + dist_max){ // discard classified voxels corresponding too far from the frontier
      for(int iter_g = 0 ; iter_g < n_groups ; iter_g++){
        while(dist_center[ls_groups[iter_g][init_groups[iter_g]]] < current_dist && (init_groups[iter_g] + 1) < size_groups[iter_g]){ // -1 to keep at least one element in each group
          init_groups[iter_g]++;
        }  
      }
      current_dist += dist_max;
    }
    
    n_neighbors = W_p[site + 1]-W_p[site] ;
    
    test_neighbor = true ; // just to launch the loop
    iter_group = -1;
    
    while(test_neighbor && iter_group < (n_groups - 1)){
      iter_group++;
      test_neighbor = false;
      
      for(int iter_neighbors = 0 ; iter_neighbors < n_neighbors; iter_neighbors++){ // iterate over sites belonging to each group to see if they are in the neighborhood of the point
        
        index_neighbor_tempo = W_i[W_p[site]+iter_neighbors];
        if(dist_center[index_neighbor_tempo] < dist_center[site]+CST_EPSILON){ // check if the neighbor point has already been classified (i.e has inferior radius)
          
          for(int iter_sitegroup = init_groups[iter_group] ; iter_sitegroup < size_groups[iter_group] ; iter_sitegroup++){
            if(abs(ls_groups[iter_group][iter_sitegroup]-index_neighbor_tempo) < 0.5){ // test whether neighborhing site have the same index than an already allocated sites
              test_neighbor = true;
            } 
          }
          
        }
      } 
    }
    
    if(test_neighbor == false){ // add site to existing group
      ls_groups[iter_group].push_back(site);
      size_groups[iter_group]++;
    }else{ // create a new group
      ls_groups.push_back(sites);
      ls_groups[n_groups].push_back(site);
      size_groups.push_back(1);
      init_groups.push_back(0);
      n_groups++;
    }
    
    if(verbose && iter_site>=value_trace){
      value_trace = min(1.0 * n, value_trace + n / CST_PROGRESS);
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
