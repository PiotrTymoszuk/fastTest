#ifndef __kappa__
#define __kappa__

#include <Rcpp.h>

using namespace Rcpp;

NumericVector kappaCpp(IntegerVector x,
                       IntegerVector y,
                       String  method);

#endif // __kappa__
