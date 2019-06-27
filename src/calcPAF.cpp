#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
using namespace Rcpp;

// [[Rcpp::export]]
void calcPAF(Rcpp::List eList, Rcpp::List alleleProbs, Rcpp::StringVector tipLabels, Rcpp::StringVector basename, Rcpp::StringVector chrom,
	     Rcpp::NumericVector start, int binSize) {
  // Set variables
  int site=1;
  int nAlleles=2;
  std::string f1=as<std::string>(basename[0])+".paf";
  std::string f2=as<std::string>(basename[0])+".bed";
  // Open file connection
  std::ofstream paffile (f1.c_str());
  std::ofstream bedfile (f2.c_str());
  // Write file header
  paffile << "site" << "\t" << "species" << "\t" << "A1" << "\t" << "A2" << std::endl;
  for(int i=0;i<eList.size();i++){
    Rcpp::LogicalVector el=eList[i];
    Rcpp::NumericMatrix aMat = alleleProbs[i];
    int nCol= aMat.ncol();
    Rcpp::NumericVector aprob(nCol);
    int eleLen = -1; 
    for(int l =0; l <el.size();l++){
      // If in an element accumulate probabilities 
      if(el(l)){
  	aprob=aprob+aMat(l,_);
	eleLen++;
      }
      // If just exiting element or at end of region and inside element, write to file
      if( (l+1 == el.size() && el(l)) || (l+1 != el.size() && el(l) && !el(l+1)) ) {
	// Write output for each species line by line
	for(int j=0; j< nCol; j=j+nAlleles){
  	  Rcpp::NumericVector a=aprob[Range(j,j+nAlleles-1)];
  	 paffile << site << "\t" << tipLabels(j/nAlleles) << "\t" << a << std::endl;
  	}
	// Write out to bed file
	int s = start(i)+((l-eleLen)*binSize);
	int e= start(i)+((l+1)*binSize);
	bedfile << chrom(i) << "\t" << s  << "\t" << e  << "\t" << site <<std::endl ;
	// Reset probability accumulator
  	Rcpp::NumericVector aprob(nCol) ;
	// Increment site count
  	site=site+1;
	// Reset element length
	eleLen = -1;
      }                     
    }            
  }
  paffile.close();
  bedfile.close();
}
