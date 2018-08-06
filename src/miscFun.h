#ifndef MISCFUN_H
#define MISCFUN_H

#include <Rcpp.h>
using namespace Rcpp;

double logSumExp(NumericVector x);
double logSumExp(double x, double y);
long double logSumExp(long double x, long double y);
double logMinusExp(double a, NumericVector x);
double logMinusExp(double x, double y);
NumericMatrix replicatePmat(NumericMatrix pmat,int nrpl);

#endif
