/*** Utilities for numeric matrices */

#include <Rcpp.h>
#include <Rmath.h>
#include "contingencyUtils.h"
#include "vectorUtils.h"
#include "transformUtils.h"

using namespace Rcpp;

// functions for retrieval of matrix diagonals, row and column
// sums, and total sums (and their frequency counterparts).
// In the sum functions, NA's are silently removed

double mtxTotalSum(NumericMatrix x) {

  // sum of all elements

  double mtx_sum = 0;
  NumericVector row_vec(x.ncol());

  for(int i = 0; i < x.nrow(); ++i) {

    row_vec = na_omit(x(i, _));

    mtx_sum += sum(row_vec);

  }

  return mtx_sum;

}

NumericVector mtxRowSum(NumericMatrix x) {

  // a numeric vector with sums of the rows

  NumericVector mtx_sum(x.nrow());
  NumericVector row_vec(x.ncol());

  for(int i = 0; i < x.nrow(); ++i) {

    row_vec = na_omit(x(i, _));

    mtx_sum[i] = sum(row_vec);

  }

  return mtx_sum;

}

NumericVector mtxRowFreq(NumericMatrix x) {

  return mtxRowSum(x)/mtxTotalSum(x);

}

NumericVector mtxColSum(NumericMatrix x) {

  // a numeric vector with sums of the rows

  NumericVector mtx_sum(x.ncol());
  NumericVector col_vec(x.nrow());

  for(int i = 0; i < x.ncol(); ++i) {

    col_vec = na_omit(x(_, i));

    mtx_sum[i] = sum(col_vec);

  }

  return mtx_sum;

}

NumericVector mtxColFreq(NumericMatrix x) {

  return mtxColSum(x)/mtxTotalSum(x);

}

NumericVector mtxDiag(NumericMatrix x) {

  // extraction of diagonal elements

  int x_ncol = x.ncol();
  int x_nrow = x.nrow();

  if(x_ncol != x_nrow) stop("The matrix is not symmetric.");

  NumericVector diag_res(x_nrow);

  for(int i = 0; i < x_ncol; ++i) {

    diag_res[i] = x(i, i);

  }

  return diag_res;

}

// absolute value, multiplication by a double value and a matrix

NumericMatrix mtxAbs(NumericMatrix x) {

  int n_row = x.nrow();

  for(int i = 0; i < n_row; ++i) {

    x(i, _) = Abs(x(i, _), false);
  }

  return x;

}

NumericMatrix mtxConstProd(NumericMatrix x, double d) {

  /// multiplies all elements of a matrix by a constant

  int n_col = x.ncol();

  for(int j = 0; j < n_col; ++j) {

    NumericVector col_vec = x(j, _);

    for(int i = 0; i < col_vec.length(); ++i) col_vec[i] *= d;

    x(j, _) = col_vec;

  }

  return x;

}

// [[Rcpp::export]]

NumericMatrix mtxProduct (NumericMatrix x, NumericMatrix y) {

  /// element-wise multiplication of two matrices of equal dimensions

  /// entry control

  if((x.ncol() != y.ncol()) | (x.nrow() != y.nrow())) {

    stop("x and y matrices must have equal dimensions");

  }

  /// result container

  int n_row = x.nrow();
  int n_col = x.ncol();

  NumericMatrix res(n_row, n_col);

  /// multiplication

  for(int i = 0; i < n_col; ++i) {

    res(_, i) = product(x(_, i), y(_, i));

  }

  return res;

}

NumericMatrix constMtxdelta (double d, NumericMatrix x) {

  int n_col = x.ncol();

  for(int j = 0; j < n_col; ++j) {

    NumericVector col_vec = x(_, j);

    for(int i = 0; i < col_vec.length(); ++i)  col_vec[i] = d - col_vec[i];

    x(_, j) = col_vec;

  }

  return x;

}

// END
