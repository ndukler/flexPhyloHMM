#include <Rcpp.h>
#include "miscFun.h"
using namespace Rcpp;

// Takes computes tree emissions for a single chain.  
// [[Rcpp::export]] 
NumericVector absorbingStateLogProb(NumericMatrix& conditionalEmisLogProb, NumericVector& dataLogProb, NumericVector& enumLogProb, double absorbLogProb ){
  NumericVector asLogProb(conditionalEmisLogProb.nrow());
  for(int i=0; i<conditionalEmisLogProb.nrow();i++){
    asLogProb(i)=logMinusExp(dataLogProb(i),conditionalEmisLogProb(i,_)+enumLogProb)-absorbLogProb;
  }
  return(asLogProb);
}
