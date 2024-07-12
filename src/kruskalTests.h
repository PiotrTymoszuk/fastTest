#ifndef __kruskalTests__
#define __kruskalTests__

#include <Rcpp.h>

using namespace Rcpp;

NumericVector kruskalVec(NumericVector x,
                         IntegerVector f,
                         bool crash);

#endif // __kruskalTests__
