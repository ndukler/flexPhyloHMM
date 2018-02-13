#include <Rcpp.h>
using namespace Rcpp;

//' Vector logSumExp
//' 
//' This function computes the sum(e^x) of a vector x without leaving log space
//'
//' @param x A numeric vector
//' @export
// [[Rcpp::export(rng = false)]]
double logSumExp(NumericVector x){
  double a=max(x);
  double out=a+log(sum(exp(x-a)));
  return(out);
}

//' Vector logMinusExp
//' 
//' This function computes the  e^a-sum(e^x) of a scalar 'a' and vector 'x' without leaving log space
//'
//' @param a A scalar
//' @param x A numeric vector
//' @export
// [[Rcpp::export(rng = false)]]
double logMinusExp(double a, NumericVector x){
  double out=a + log(std::max(0.0,1 - sum(exp(x - a))));
  return(out);
}

NumericMatrix replicatePmat(NumericMatrix pmat,int nrpl){
  if(pmat.ncol()==1){
    NumericMatrix tmp(pmat.nrow(),nrpl);
    for(int i=0;i<nrpl;i++){
      tmp(_,i)=pmat(_,0);
    }
    pmat=Rcpp::clone(tmp);
  }
  return(pmat);
}
