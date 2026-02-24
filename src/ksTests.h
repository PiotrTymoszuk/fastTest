#ifndef __ksTests__
#define __ksTests__

#include <Rcpp.h>

using namespace Rcpp;

NumericVector ksTestCpp(NumericVector x,
                        NumericVector y,
                        String alternative);

NumericVector ksTestVec(NumericVector x,
                        IntegerVector f,
                        String alternative);


#endif // __ksTests__
