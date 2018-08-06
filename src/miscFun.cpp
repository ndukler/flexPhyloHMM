#include <Rcpp.h>
#include <cmath> 
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

//' Vector logSumExpHP
//' 
//' This function computes the sum(e^x) of a vector x without leaving log space in high precision using an iterative calculation
//'
//' @param x A numeric vector
//' @export
// [[Rcpp::export(rng = false)]]
double logSumExpHP(NumericVector z){
  double mx = max(z);
  NumericVector x = exp(z-mx);
  double ks =0;
  double c=0;
  for(int i =0;i<x.size();i++){
    double y = x(i)-c;
    double kt = ks+y;
    c=(kt-ks)-y;
    ks=kt;
  }
  double out = mx + log(ks);
  return(out);
}

double logSumExp(double x, double y){
  double a = std::max(x,y);
  double b = std::min(x,y);
  double out=a+log1p(std::exp(b-a));
  return(out);
}

long double logSumExp(long double x, long double y){
  long double a = std::max(x,y);
  long double b = std::min(x,y);
  long double out=a+log1p(std::exp(b-a));
  return(out);
}

double logMinusExp(double x, double y){
  double out;
  if(y>x){
    out=R_NaN;
  } else { 
    out=x + log1p(-exp(y-x));
  } 
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
double logMinusExp(double x, NumericVector y){
  double out;
  double interSum = logSumExp(y);
  if(interSum>x){
    out=R_NaN;
    double ovr = logMinusExp(interSum,x);
    Rcpp::Rcout.precision(10);
    Rcpp::Rcout << "X: " << x << "  Y: " << interSum << "  Over-run:"  << ovr << "  Over-Run - X: " << ovr-x  <<std::endl;
    // throw std::range_error("Inadmissible value");
  } else { 
    out=x + log1p(-exp(interSum-x));
  } 
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
