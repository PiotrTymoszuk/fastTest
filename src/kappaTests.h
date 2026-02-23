#ifndef __kappaTests__
#define __kappaTests__

#include <Rcpp.h>

using namespace Rcpp;

// permutation tests

NumericVector permKappaVec(IntegerVector x,
                           IntegerVector y,
                           String method,
                           String alternative,
                           int n_iter);

// bootstrap tests

NumericVector bootKappaVec(IntegerVector x,
                           IntegerVector y,
                           String method,
                           String ci_type,
                           double conf_level,
                           int n_iter);

#endif // __kappaTests__
