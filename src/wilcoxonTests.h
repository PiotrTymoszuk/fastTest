#ifndef __wilcoxonTests__
#define __wilcoxonTests__

#include <Rcpp.h>

using namespace Rcpp;

// effect size

NumericVector biserialR(NumericVector x,
                        IntegerVector f,
                        bool paired,
                        bool hl_estimate);

// Wilcoxon tests

NumericVector testWilcoxonVec(NumericVector x,
                              IntegerVector f,
                              String type,
                              String alternative,
                              bool exact,
                              bool correct,
                              bool hl_estimate,
                              bool crash);


#endif // __wilcoxonTests__
