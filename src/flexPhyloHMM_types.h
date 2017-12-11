#ifndef DATTYPES_H
#define DATTYPES_H

#include <Rcpp.h>

// Define type to pass list of numeric matricies
typedef Rcpp::ListOf<Rcpp::NumericMatrix> NumMatList;
// Define type to pass list of numeric vectors
typedef Rcpp::ListOf<Rcpp::NumericVector> NumVecList;

#endif
