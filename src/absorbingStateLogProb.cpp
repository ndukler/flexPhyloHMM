#include <Rcpp.h>
#include "miscFun.h"
using namespace Rcpp;

// Takes computes tree emissions for a single chain.  
// [[Rcpp::export]] 
NumericVector absorbingStateLogProb(NumericMatrix& conditionalEmisLogProb, NumericVector& dataLogProb, NumericVector& enumLogProb, double absorbLogProb ){
  NumericVector asLogProb(conditionalEmisLogProb.nrow());
  for(int i=0; i<conditionalEmisLogProb.nrow();i++){
    double jitter = logSumExp(dataLogProb(i),dataLogProb(i)-27.0);
    // double jitter = dataLogProb(i);
    asLogProb(i)=logMinusExp(jitter,conditionalEmisLogProb(i,_)+enumLogProb)-absorbLogProb;
  }
  return(asLogProb);
}
