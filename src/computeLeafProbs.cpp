#include <Rcpp.h>
#include "flexPhyloHMM_types.h"
#include "miscFun.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix computeLeafProb(NumMatList& data, NumMatList& scaleFactors, NumMatList& deletionRanges, 
					   NumMatList muBack, NumMatList mixPeak, NumVecList weights, NumVecList sampleNormFactor, NumericVector dispParams) {
  int nspecies=data.size();
  int nobs=data[0].nrow();
  int ncomp=weights[0].size();
  // Initialize matrix to hold 0/1 state probabilities for each species, two columns per species, in the same order as the data list
  NumericMatrix out(nobs,nspecies*2);
  // Iterate over species
  for(int i=0;i<nspecies;i++){
    int nrpl = data[i].ncol();
    int delInd=0;
    int sCol=i*2;
    // Expand parameter matricies if more replicates than columns in the parameter matrix
    NumericMatrix spMuBack=replicatePmat(muBack[i],nrpl);
    NumericMatrix spMixPeak=replicatePmat(mixPeak[i],nrpl);
    // Scale means by sample normalization factors
    for(int r=0;r<nrpl;r++){
      spMixPeak(_,r)=spMixPeak(_,r)*sampleNormFactor[i](r);
      spMuBack(_,r)=spMuBack(_,r)*sampleNormFactor[i](r);
    }
    // First compute the log-probabilities under both negative binomial models

    // Iterate over observations
    for(int j=0;j<nobs;j++){
      if(j >= deletionRanges[i](delInd,0) && j < deletionRanges[i](delInd,1)){ // if in a deletion region just set probability of a peak to -Inf  and log probability of no element to 0
	out(j,sCol)=0;
	out(j,sCol+1)=- std::numeric_limits<double>::infinity();
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
	    out(j,sCol) = out(j,sCol) + R::dnbinom_mu(data[i](j,k),expSizeBack,expMuBack,true);
	    // Compute the probability for the active state, summing over all components	
	    for(int m=0;m<ncomp;m++){
	      double expMu=spMixPeak(m,k)*scaleFactors[i](j,0);
	      double expSize= 1.0/(dispParams(0)+dispParams(1)/expMu);
	      cProb(m) = cProb(m) + R::dnbinom_mu(data[i](j,k),expSize,expMu,true);
	    }
	  } // exiting iteration over replicates
	  out(j,sCol+1)=logSumExp(cProb+log(weights[i]));
	}
	else { // If the scale factor is zero just set probability of active state to 0.01 (1/10th of inactive state)
	  out(j,sCol)=0;
	  out(j,sCol+1)=0;
	}
      }
    }
  }
  return(out);
}
