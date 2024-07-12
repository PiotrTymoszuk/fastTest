#ifndef __oneBlockAnova__
#define __oneBlockAnova__

#include <Rcpp.h>

using namespace Rcpp;

NumericVector oneBlockAnovaVec(NumericVector x,
                               IntegerVector f,
                               IntegerVector b);


#endif // __oneBlockAnova__
