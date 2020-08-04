#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix computeProb(double lambdamin,double lambdamax,double beta,NumericMatrix matpval,Function integral) {
//Rcpp::Function integral=Environment::global_env()["integral"];
  int nbrow= matpval.nrow();
  int nbcol= matpval.ncol();
  NumericMatrix res (nbrow,nbcol);
  double pval;
 // beta= Environment::global_env()["beta"];
 //lambdamin=Environment::global_env()["lambdamin"];
  //lambdamax=Environment::global_env()["lambdamax"];
  for (int i=0;i<nbrow;i++)
  {
   
    for (int j=0;j<nbcol;j++)
    {
     pval=matpval(i,j);
     
      if (i==j)continue;
      res(i,j)=1/Rcpp::as<double>(integral(beta,pval,lambdamin,lambdamax));
      //{
      //  res(i,j)=0;
     //}
     // else
      //{
        
        
        
       // 
      //}
     
    }
  }
  return (res);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


