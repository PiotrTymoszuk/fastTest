#ifndef __rankUtils__
#define __rankUtils__

#include <Rcpp.h>

using namespace Rcpp;

// rank calculation

NumericVector Rank(NumericVector x,
                   String ties_method,
                   bool na_rm);

NumericVector signedRank(NumericVector x,
                         String ties_method,
                         bool na_rm);

// ranking of vectors

List rankVectors(NumericVector x, IntegerVector f);

List rankBlockVectors(NumericVector x,
                      IntegerVector f,
                      IntegerVector b);

#endif // __rankUtils__
