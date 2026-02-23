#ifndef __vectorUtils__
#define __vectorUtils__

#include <Rmath.h>

using namespace Rcpp;

// functions for operations on vectors

NumericVector delta(NumericVector x, NumericVector y);
NumericVector product(NumericVector x, NumericVector y);

NumericVector Sign(NumericVector x, bool na_rm);
NumericVector Abs(NumericVector x, bool na_rm);

// outer and cross operations

NumericMatrix outerSum(NumericVector x, NumericVector y, bool na_rm);
NumericMatrix outerDelta(NumericVector x, NumericVector y, bool na_rm);
NumericMatrix outerProduct(NumericVector x, NumericVector y, bool na_rm);

NumericMatrix crossProduct(NumericVector x, NumericVector y);

// vectors of values of a matrix

NumericVector matrix2vector(NumericMatrix x, bool by_row);

// sorting and ordering, duplicates

NumericVector orderVector(NumericVector x, NumericVector y);
IntegerVector countMulti(NumericVector x);

// resampling

NumericMatrix resampleVec(NumericVector x, int n_iter, bool replace);

// string and numeric sequences

CharacterVector stringSeq(String prefix, int start, int end);
NumericVector intSeq(int start, int end);

#endif // __vectorUtils__
