#ifndef __transformUtils__
#define __transformUtils__

#include <Rcpp.h>

using namespace Rcpp;

// checking if R object have names

bool checkListNames(const List &x);
LogicalVector checkNames (const NumericMatrix &x);

// exposed functions for transformation of vectors and matrices

List completeCases(NumericVector x, IntegerVector f);

List Split(NumericVector x, IntegerVector f);
List SplitMtx(NumericMatrix x, IntegerVector f);

NumericVector Concat(NumericVector x, NumericVector y);

// quality control for the input vectors in T and Wilcoxon tests

List procVectors(NumericVector x,
                 NumericVector y,
                 bool paired);

// checking sizes of categories of a splittin factor

bool checkEqualSplits(NumericVector x);

#endif // __transformUtils__
