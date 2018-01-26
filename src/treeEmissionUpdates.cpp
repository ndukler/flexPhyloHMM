#include <Rcpp.h>
#include "flexPhyloHMM_types.h"
#include "miscFun.h"
using namespace Rcpp;

// Takes computes tree emissions for a single chain.  
// [[Rcpp::export]] 
NumericMatrix updateTreeEmisProbCpp(NumMatList& data, NumMatList& scaleFactors, NumericMatrix treeStates, int nstates, NumMatList& mu, NumMatList& size) {
  int nspecies=data.size();
  int nobs=data[0].nrow();
  // Match parameterizations
  // NumericVector prob = size/(size+mu);
  NumericMatrix out(nobs,nstates);
  // Iterate over species
  for(int i=0;i<nspecies;i++){
    int nrpl = data[i].ncol();
    // First compute the log-probabilities under both negative binomial models
    NumericMatrix logProb(nobs,2);
    // Iterate over observations
    for(int j=0;j<nobs;j++){
      // Sum over replicates
      for(int k=0;k<nrpl;k++){
  	logProb(j,0) = logProb(j,0) + R::dnbinom_mu(data[i](j,k),size[i](0,k),mu[i](0,k)*scaleFactors[i](j,0),true);
  	logProb(j,1) = logProb(j,1) + R::dnbinom_mu(data[i](j,k),size[i](1,k),mu[i](1,k)*scaleFactors[i](j,0),true);
      }
    }
    // Now iterate over states
    for(int l=0;l<nstates;l++){
      for(int j=0;j<nobs;j++){
	out(j,l)=out(j,l)+logProb(j,treeStates(l,i));
      }
    }
  }
  return(out);
}


// Takes computes tree emissions for a single chain.  
// [[Rcpp::export]] 
NumericMatrix updateTreeEmisProbDTCpp(NumMatList& data, NumMatList& scaleFactors, NumericMatrix treeStates, int nstates, NumMatList mu, NumericVector dispParams) {
  int nspecies=data.size();
  int nobs=data[0].nrow();
  // Match parameterizations
  // NumericVector prob = size/(size+mu);
  NumericMatrix out(nobs,nstates);
  // Iterate over species
  for(int i=0;i<nspecies;i++){
    int nrpl = data[i].ncol();
    // First compute the log-probabilities under both negative binomial models
    NumericMatrix logProb(nobs,2);
    // Iterate over observations
    for(int j=0;j<nobs;j++){
      // Sum over replicates
      for(int k=0;k<nrpl;k++){
	double expMuBack=mu[i](0,k)*scaleFactors[i](j,0);
	double expMuPeak=mu[i](1,k)*scaleFactors[i](j,0);
	double expSizeBack= 1.0/(dispParams(0)+dispParams(1)/expMuBack);
	double expSizePeak= 1.0/(dispParams(0)+dispParams(1)/expMuPeak);
  	logProb(j,0) = logProb(j,0) + R::dnbinom_mu(data[i](j,k),expSizeBack,expMuBack,true);
  	logProb(j,1) = logProb(j,1) + R::dnbinom_mu(data[i](j,k),expSizePeak,expMuPeak,true);
      }
    }
    // Now iterate over states
    for(int l=0;l<nstates;l++){
      for(int j=0;j<nobs;j++){
	out(j,l)=out(j,l)+logProb(j,treeStates(l,i));
      }
    }
  }
  return(out);
}

// Takes computes tree emissions for a single chain.  
// [[Rcpp::export]] 
NumericMatrix updateTreeEmisProbDTWindowCpp(NumMatList& data, NumMatList& scaleFactors, NumericMatrix treeStates, int nstates, NumMatList mu, NumericVector dispParams, int winSize) {
  int nspecies=data.size();
  int nobs=data[0].nrow();
  // Match parameterizations
  // NumericVector prob = size/(size+mu);
  NumericMatrix out(nobs,nstates);
  // Iterate over species
  for(int i=0;i<nspecies;i++){
    int nrpl = data[i].ncol();
    // First compute the log-probabilities under both negative binomial models
    NumericMatrix logProb(nobs,2);
    // Iterate over observations
    for(int j=0;j<nobs;j++){
      // Sum over replicates
      for(int k=0;k<nrpl;k++){
	double expMuBack=mu[i](0,k)*scaleFactors[i](j,0);
	double expMuPeak=mu[i](1,k)*scaleFactors[i](j,0);
	double expSizeBack= 1.0/(dispParams(0)+dispParams(1)/expMuBack);
	double expSizePeak= 1.0/(dispParams(0)+dispParams(1)/expMuPeak);
  	logProb(j,0) = logProb(j,0) + R::dnbinom_mu(data[i](j,k),expSizeBack,expMuBack,true);
  	logProb(j,1) = logProb(j,1) + R::dnbinom_mu(data[i](j,k),expSizePeak,expMuPeak,true);
      }
    }
    // EXPERIMENTAL: Compute Log prob in windows
    //  Set window size
    // Now iterate over states
    for(int l=0;l<nstates;l++){
      for(int j=0;j<nobs;j++){
	//sum ober window
	for(int w=std::max(0,j-winSize); w < std::min(nobs,j+winSize);w++){
	  out(j,l)=out(j,l)+logProb(w,treeStates(l,i));
	}
      }
    }
  }
  return(out);
}


// Takes computes tree emissions for a single chain under a gamma mixure model for the peak state
// [[Rcpp::export]]
NumericMatrix updateTreeEmisProbGMixCpp(NumMatList& data, NumMatList& scaleFactors, NumericMatrix treeStates, int nstates, NumMatList muBack, NumMatList mixPeak, NumericVector weights, NumericVector dispParams) {
  int nspecies=data.size();
  int nobs=data[0].nrow();
  int ncomp=weights.size();
  // Initialize logProb output matrix
  NumericMatrix out(nobs,nstates);
  // Iterate over species
  for(int i=0;i<nspecies;i++){
    int nrpl = data[i].ncol();
    // First compute the log-probabilities under both negative binomial models
    NumericMatrix logProb(nobs,2);
    // Iterate over observations
    for(int j=0;j<nobs;j++){
      // Sum over replicates
      for(int k=0;k<nrpl;k++){
	// Compute the probability for the background state
	double expMuBack=muBack[i](0,k)*scaleFactors[i](j,0);
	double expSizeBack= 1.0/(dispParams(0)+dispParams(1)/expMuBack);
	logProb(j,0) = logProb(j,0) + R::dnbinom_mu(data[i](j,k),expSizeBack,expMuBack,true);
	// Compute the probability for the active state, summing over all components	
	double mixDen=0;
	for(int m=0;m<ncomp;m++){
	  double expMu=mixPeak[i](m,k)*scaleFactors[i](j,0);
	  double expSize= 1.0/(dispParams(0)+dispParams(1)/expMu);
	  mixDen = mixDen  + weights(m) * R::dnbinom_mu(data[i](j,k),expSize,expMu,false);
	}
	logProb(j,1)+=std::log(mixDen);
      }
    }
    // Now iterate over states
    for(int l=0;l<nstates;l++){
      for(int j=0;j<nobs;j++){
	out(j,l)=out(j,l)+logProb(j,treeStates(l,i));
      }
    }
  }
  return(out);
}



// [[Rcpp::export]]
NumericMatrix updateTreeEmisProbFlexMixCpp(NumMatList& data, NumMatList& scaleFactors, NumMatList& deletionRanges,NumericMatrix treeStates, int nstates, 
					   NumMatList muBack, NumMatList mixPeak, NumVecList weights, NumVecList sampleNormFactor, NumericVector dispParams) {
  int nspecies=data.size();
  int nobs=data[0].nrow();
  int ncomp=weights[0].size();
  // Initialize logProb output matrix
  NumericMatrix out(nobs,nstates);
  // Iterate over species
  for(int i=0;i<nspecies;i++){
    int nrpl = data[i].ncol();
    int delInd=0;
    // Expand parameter matricies if more replicates than columns in the parameter matrix
    NumericMatrix spMuBack=replicatePmat(muBack[i],nrpl);
    NumericMatrix spMixPeak=replicatePmat(mixPeak[i],nrpl);
    // Scale means by sample normalization factors
    for(int r=0;r<nrpl;r++){
      spMixPeak(_,r)=spMixPeak(_,r)*sampleNormFactor[i](r);
      spMuBack(_,r)=spMuBack(_,r)*sampleNormFactor[i](r);
    }
    // First compute the log-probabilities under both negative binomial models
    NumericMatrix logProb(nobs,2);
    // Iterate over observations
    for(int j=0;j<nobs;j++){
      if(j >= deletionRanges[i](delInd,0) && j < deletionRanges[i](delInd,1)){ // if in a deletion region just set probability of a peak to -Inf  and log probability of no element to 0
	logProb(j,0)=0;
	logProb(j,1)=- std::numeric_limits<double>::infinity();
      }
      else { // If not in a deletion do normal calculation
	if(j == deletionRanges[i](delInd,1)){
	  delInd=delInd+1;
	}	
	if(scaleFactors[i](j,0)!=0){ // if scale factor is not zero do normal calculation
	  NumericVector cProb(ncomp);
	  // Sum over replicates
	  for(int k=0;k<nrpl;k++){
	    // Compute the probability for the background state
	    double expMuBack=spMuBack(0,k)*scaleFactors[i](j,0);
	    double expSizeBack= 1.0/(dispParams(0)+dispParams(1)/expMuBack);
	    logProb(j,0) = logProb(j,0) + R::dnbinom_mu(data[i](j,k),expSizeBack,expMuBack,true);
	    // Compute the probability for the active state, summing over all components	
	    for(int m=0;m<ncomp;m++){
	      double expMu=spMixPeak(m,k)*scaleFactors[i](j,0);
	    double expSize= 1.0/(dispParams(0)+dispParams(1)/expMu);
	    cProb(m) = cProb(m) + R::dnbinom_mu(data[i](j,k),expSize,expMu,true);
	    //if(isnan(R::dnbinom_mu(data[i](j,k),expSize,expMu,true))){
	    //  Rcout << data[i](j,k) << m << " " <<  k << " "  << expMu << " " << spMixPeak(m,k) << " " << scaleFactors[i](j,0) << std::endl;
	    // }
	    }
	  } // exiting iteration over replicates
	  logProb(j,1)=logSumExp(cProb+log(weights[i]));
	}
	else { // If the scale factor is zero just set probability of active state to 0.01 (1/10th of inactive state)
	  logProb(j,1)=0;
	  logProb(j,0)=0;
	}
      }
    }
    // exiting calculation of probability of peak/no-peak for species i over all positions
    // Now iterate over states
    for(int l=0;l<nstates;l++){
      for(int j=0;j<nobs;j++){
	out(j,l)=out(j,l)+logProb(j,treeStates(l,i));
      }
    }
  }
  return(out);
}
