#include <Rcpp.h>
using namespace Rcpp;

double logSumExp(NumericVector x){
  double a=max(x);
  double out=a+log(sum(exp(x-a)));
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
