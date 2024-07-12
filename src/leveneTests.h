#ifndef __leveneTests__
#define __leveneTests__

#include <Rcpp.h>

using namespace Rcpp;

// base version of Levene test for lists of numeric vectors

NumericVector leveneBase(List x, String type);

// version of the test for a numeric vector and a splitting factor

NumericVector leveneVec(NumericVector x,
                        IntegerVector f,
                        String type,
                        bool crash);

#endif // __leveneTests__
