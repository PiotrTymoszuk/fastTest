#ifndef __contingencyUtils__
#define __contingencyUtils__

#include <Rcpp.h>

using namespace Rcpp;

// exposed contingency table utilities

IntegerVector Table(NumericVector x);
NumericMatrix xTable(NumericVector x, IntegerVector f);
NumericMatrix irrTable(IntegerVector x, IntegerVector y);

#endif // __contingencyUtils__
