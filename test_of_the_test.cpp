#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// [[Rcpp::export]]
NumericVector timesTwo(NumericMatrix a, NumericMatrix b) {
  NumericMatrix aaa(a.ncol(), b.nrow());
  for(int i =0; i<a.ncol(); i++){
    aaa(j,i) = a(j,i) * b(i,j);
  }
  return(aaa);

}
