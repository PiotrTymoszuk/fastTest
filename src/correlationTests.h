#ifndef __correlationTests__
#define __correlationTests__

#include <Rcpp.h>

using namespace Rcpp;

// permutation tests

NumericVector permCorVec(NumericVector x,
                         NumericVector y,
                         String method,
                         String alternative,
                         int n_iter);

NumericVector bootCorVec(NumericVector x,
                         NumericVector y,
                         String method,
                         String ci_type,
                         double conf_level,
                         int n_iter);

#endif // __correlationTests__
