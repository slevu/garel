// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;//

//[[Rcpp::export]]
arma::vec d2root_cpp(const int n, const arma::ivec node2parent, const arma::vec node2el){
  int i, u, v;
  vec d2rs = zeros<arma::vec>(n);
  for (i = 0; i < n; i++){
    u = node2parent(i)-1;
    v = i;
    while( u >= 0){
      d2rs(i) += node2el(v);
      v = u;
      u = node2parent(v)-1;
    }
  }
  return d2rs;
}
