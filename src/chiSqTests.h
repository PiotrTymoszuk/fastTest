#ifndef __chiSqTests__
#define __chiSqTests__

#include <Rcpp.h>

using namespace Rcpp;

bool ctgFlaws(NumericMatrix ctg);

NumericVector chiSqTestTbl(NumericMatrix ctg,
                           bool correct,
                           bool crash);

NumericVector chiSqTestVec(NumericVector x,
                           IntegerVector f,
                           bool correct,
                           bool crash);

#endif // __chiSqTests__
