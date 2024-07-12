#ifndef __friedmanTests__
#define __friedmanTests__

#include <Rcpp.h>

using namespace Rcpp;

NumericVector friedmanVec(NumericVector x,
                          IntegerVector f,
                          IntegerVector b,
                          bool crash);

#endif // __friedmanTests__
