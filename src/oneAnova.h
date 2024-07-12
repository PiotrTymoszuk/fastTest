#ifndef __oneAnova__
#define __oneAnova__

#include <Rcpp.h>

using namespace Rcpp;

NumericVector oneAnovaBase(List x);

NumericVector oneAnovaVec(NumericVector x,
                          IntegerVector f,
                          bool crash);

#endif // __oneAnova__
