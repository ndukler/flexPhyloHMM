#ifndef MISCFUN_H
#define MISCFUN_H

#include <Rcpp.h>
using namespace Rcpp;

double logSumExp(NumericVector x);
NumericMatrix replicatePmat(NumericMatrix pmat,int nrpl);

#endif
