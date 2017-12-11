#include <Rcpp.h>
#include "miscFun.h"
using namespace Rcpp;

// Takes computes peak emissions for a single chain. 
 // [[Rcpp::export]] 
NumericMatrix updateNBEmisProbCpp(NumericMatrix& data, NumericVector& scaleFactors, int& nstates, NumericMatrix& muBack,NumericMatrix& muPeak, NumericMatrix& size) {
  int nobs=data.nrow();
  int nrpl=data.ncol();
  // If parameters are not not different between replicates, match column dimensionality
  muPeak=replicatePmat(muPeak,nrpl);
  muBack=replicatePmat(muBack,nrpl);
  size=replicatePmat(size,nrpl);
  NumericMatrix out(nobs,nstates);
  // Iterate over states
  for(int i=0;i<nstates;i++){
    NumericMatrix *mu;
    if(i==0){
      mu = &muBack;
    } 
    else {
      mu = &muPeak;
    }
    // Iterate over observations
    for(int j=0;j<nobs;j++){
      // Iterate over replicates
      for(int k=0;k<nrpl;k++){
  	out(j,i) = out(j,i) + R::dnbinom_mu(data(j,k),size(i,k),(*mu)(0,k)*scaleFactors(j),true);
      }
    }
  }
  return(out);
}

// Takes computes peak emissions for a single chain. 
 // [[Rcpp::export]] 
NumericMatrix updateNBEmisProbDTCpp(NumericMatrix& data, NumericVector& scaleFactors, int& nstates, NumericMatrix& muBack,NumericMatrix& muPeak, NumericVector& dispPar) {
  int nobs=data.nrow();
  int nrpl=data.ncol();
  // If parameters are not not different between replicates, match column dimensionality
  muPeak=replicatePmat(muPeak,nrpl);
  muBack=replicatePmat(muBack,nrpl);
  NumericMatrix out(nobs,nstates);
  // Iterate over states
  for(int i=0;i<nstates;i++){
    NumericMatrix *mu;
    if(i==0){
      mu = &muBack;
    } 
    else {
      mu = &muPeak;
    }
    // Iterate over observations
    for(int j=0;j<nobs;j++){
      // Iterate over replicates
      for(int k=0;k<nrpl;k++){
	double expMu=(*mu)(0,k)*scaleFactors(j);
	double expSize= 1.0/(dispPar(0)+dispPar(1)/expMu);
  	out(j,i) = out(j,i) + R::dnbinom_mu(data(j,k),expSize,expMu,true);
      }
    }
  }
  return(out);
}

// Takes computes peak emissions for a single chain. 
 // [[Rcpp::export]] 
NumericMatrix updatePeakFlexMixCpp(NumericMatrix& data, NumericVector& scaleFactors, int& nstates, NumericMatrix& muBack,NumericMatrix& mixPeak, NumericVector weights,
				   NumericVector dispPar,  NumericVector sampleNormFactor) {
  int nobs=data.nrow();
  int nrpl=data.ncol();
  int ncomp=mixPeak.nrow();
  // If parameters are not not different between replicates, match column dimensionality
  mixPeak=replicatePmat(mixPeak,nrpl);
  muBack=replicatePmat(muBack,nrpl);
  // sizeBack=replicatePmat(sizeBack,nrpl);
  // Scale means by sample normalization factors
  for(int i=0;i<nrpl;i++){
    mixPeak(_,i)=mixPeak(_,i)*sampleNormFactor(i);
    muBack(_,i)=muBack(_,i)*sampleNormFactor(i);
  }
  // Multiply the peak and back parameter by the sample scale factor
  NumericMatrix out(nobs,nstates);
  // Iterate over states
  for(int i=0;i<nstates;i++){
    // Iterate over observations
    for(int j=0;j<nobs;j++){
      // Iterate over replicates
      NumericVector cProb(ncomp);
	// if inactive, just compute probability of emission under single mu
      if(i==0){
	for(int k=0;k<nrpl;k++){
	  double expMu=muBack(0,k)*scaleFactors(j);
	  double expSize= 1.0/(dispPar(0)+dispPar(1)/expMu);
	  out(j,i) = out(j,i) + R::dnbinom_mu(data(j,k),expSize,expMu,true);
	}
      }
      else { // if active, just compute probability of emission under mixture of mu's.
	// Iterate over componenets
	for(int k=0;k<nrpl;k++){
	  for(int m=0;m<ncomp;m++){
	    double expMu=mixPeak(m,k)*scaleFactors(j);
	    double expSize= 1.0/(dispPar(0)+dispPar(1)/expMu);
	    cProb(m) = cProb(m) + R::dnbinom_mu(data(j,k),expSize,expMu,true);
	  }	  
	}
	out(j,i)=logSumExp(cProb+log(weights));
      }
    }
  }
  return(out);
}

