#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix matmult(NumericMatrix x, NumericMatrix y) {
  NumericMatrix ans(x.nrow(), y.ncol());
  
  for(int i = 0; i<x.nrow(); i++){
    for(int j=0; j<y.ncol(); j++){
      for(int k =0; k<y.nrow(); k++){
        
        ans(i,j) += x(i,k) * y(k,j);
      }
    }}
  
  return ans;
}
