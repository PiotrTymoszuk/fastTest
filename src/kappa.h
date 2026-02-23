#ifndef __kappa__
#define __kappa__

#include <Rcpp.h>

using namespace Rcpp;

NumericVector kappaCpp(IntegerVector x,
                       IntegerVector y,
                       String type);

#endif // __kappa__
