#include <Rcpp.h>
#include "flexPhyloHMM_types.h"
#include "miscFun.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector treeLL(const NumericMatrix& leafProb, const NumericMatrix& qConcat, const NumericMatrix& traversal, const NumericVector& tips, const NumericVector& logBaseFreq) {
  int sites= leafProb.nrow();
  NumericVector asLogProb(sites); // Numeric vector of the for the probability of the absorbing state
  int nNode=Rcpp::max(traversal(_,0))+1; // The total number of nodes on the tree
  int root=traversal(traversal.nrow()-1,0); // Get the index of the root node
  // loop over sites
  for(int i=0;i<sites;i++){
    NumericMatrix nodeLogProb(nNode,2);
    // Normally here's where the 0/1 probability calculations from the data go for each tip of the tree
    for(int n=0;n<tips.size();n++){
      nodeLogProb(tips(n),0)=leafProb(i,tips(n)*2);
      nodeLogProb(tips(n),1)=leafProb(i,tips(n)*2+1);
    }
    // Now compute the probability of 0/1 for the interior nodes
    for(int n=0;n<traversal.nrow();n++){
      int parentInd=traversal(n,0);
      int childInd=traversal(n,1);
      int qRow=traversal(n,2)*2; // The starting row for the edge
      
      nodeLogProb(parentInd,0) = nodeLogProb(parentInd,0)+ logSumExp(NumericVector::create(qConcat(qRow,0) + nodeLogProb(childInd,0),qConcat(qRow,1) + nodeLogProb(childInd,1)));
      nodeLogProb(parentInd,1) = nodeLogProb(parentInd,1) + logSumExp(NumericVector::create(qConcat(qRow+1,0) + nodeLogProb(childInd,0),qConcat(qRow+1,1) + nodeLogProb(childInd,1)));

    }
    asLogProb(i)=logSumExp(nodeLogProb(root,_)+logBaseFreq);
  }
  return(asLogProb);
}
