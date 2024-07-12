#ifndef __numericUtils__
#define __numericUtils__

#include <Rcpp.h>

using namespace Rcpp;

// effect size

NumericVector cohenD(NumericVector x,
                     NumericVector y,
                     String type);

// T tests

NumericVector tTestStd(NumericVector x,
                       NumericVector y,
                       String type,
                       String alternative,
                       double conf_level,
                       bool crash);

NumericVector tTestVec(NumericVector x,
                       IntegerVector f,
                       String type,
                       String alternative,
                       double conf_level,
                       bool crash);


#endif // __numericUtils__
