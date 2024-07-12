#ifndef __contingencyUtils__
#define __contingencyUtils__

#include <Rcpp.h>

using namespace Rcpp;

// exposed contingency table utilities

IntegerVector Table(NumericVector x);
NumericMatrix xTable(NumericVector x, IntegerVector f);

#endif // __contingencyUtils__
