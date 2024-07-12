#ifndef __metaEstimates__
#define __metaEstimates__

#include <Rcpp.h>

using namespace Rcpp;

// fixed and random effect meta-estimates: helper functions

NumericVector MetaFun(NumericVector y,
                      NumericVector e,
                      double tausq,
                      String alternative,
                      double conf_level);

NumericVector estimateVar(NumericVector y, NumericVector e);

// computation of meta estimates for various inputs

NumericVector metaVec(NumericVector y,
                      NumericVector e,
                      String type,
                      String alternative,
                      double conf_level,
                      bool crash);


#endif // __metaEstimates__
