#ifndef __matrixUtils__
#define __matrixUtils__

#include <Rcpp.h>

using namespace Rcpp;

// exposed matrix utilities

/// frequencies, sums, and diagnoal elements

double mtxTotalSum(NumericMatrix x);

NumericVector mtxRowSum(NumericMatrix x);
NumericVector mtxRowFreq(NumericMatrix x);

NumericVector mtxColSum(NumericMatrix x);
NumericVector mtxColFreq(NumericMatrix x);

NumericVector mtxDiag(NumericMatrix x);

/// algebra

NumericMatrix mtxAbs(NumericMatrix x);
NumericMatrix mtxConstProd(NumericMatrix x, double d);
NumericMatrix mtxProduct (NumericMatrix x, NumericMatrix y);
NumericMatrix constMtxdelta (double d, NumericMatrix x);

#endif // __matrixUtils__
