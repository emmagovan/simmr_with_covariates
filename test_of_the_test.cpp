#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// [[Rcpp::export]]
NumericVector timesTwo(NumericMatrix expf, int n_sources, int n) {
  NumericVector sumexpf(n_sources);
  
  for(int i = 0; i<n; i++){
    for(int j=0; j<n_sources; j++){
    sumexpf(i) +=expf(i,j);
  }
  }
  return(sumexpf);
}
