/*** Cohen's kappa inter-rate reliability statistic*/

#include <Rcpp.h>
#include <Rmath.h>
#include "contingencyUtils.h"
#include "vectorUtils.h"
#include "transformUtils.h"
#include "matrixUtils.h"

using namespace Rcpp;

// Cohen's kappa for a pair of integer vectors

// [[Rcpp::export]]

NumericVector kappaCpp(IntegerVector x,
                       IntegerVector y,
                       String method) {

  /// contingency matrix: length checks and NA removal
  /// is done by the internal utility

  NumericMatrix ctg_mtx = irrTable(x, y);

  int n_col = ctg_mtx.ncol();

  // the output vector and observation number checks

  NumericVector res(2, NA_REAL);
  res.names() = CharacterVector({"n", "kappa"});

  double ctg_total_freq = mtxTotalSum(ctg_mtx);

  res[0] = ctg_total_freq;

  if(n_col < 2) {

    warning("At least two non-NA observations are required");

    return res;

  }

  /// matrix diagonal and margins

  NumericVector ctg_diag = mtxDiag(ctg_mtx);

  NumericVector ctg_row_freq = mtxRowFreq(ctg_mtx);
  NumericVector ctg_col_freq = mtxColFreq(ctg_mtx);

  /// po: observed agreement (matrix diagonal by total frequency)
  /// pc: theoretical coincidence by chance
  /// used for unweighted kappa

  double po;
  double pc;

  if(method == "unweighted") {

    po = sum(ctg_diag)/ctg_total_freq;
    pc = crossProduct(ctg_col_freq, ctg_row_freq)(0, 0);

  } else {

    /// common variables for calculation of the weights

    NumericVector dim_seq = intSeq(1, n_col);

    NumericMatrix outer_delta = mtxAbs(outerDelta(dim_seq, dim_seq, false));

    outer_delta = mtxConstProd(outer_delta, 1.0/(n_col - 1));

    NumericMatrix weights;

    if(method == "equal") {

      /// equal-spaced weights

      weights = constMtxdelta(1.0, outer_delta);

    } else {

      /// Fleiss-Cohen weights

      weights = constMtxdelta(1.0,
                              mtxProduct(outer_delta,
                                         outer_delta));

    }

    po = mtxTotalSum(mtxProduct(weights, ctg_mtx))/ctg_total_freq;

    pc = mtxTotalSum(mtxProduct(weights,
                                outerProduct(ctg_col_freq,
                                             ctg_row_freq,
                                             false)));

  }

  /// the output vector

  res[1] = (po - pc)/(1 - pc);

  return res;

}

// Cohen's kappa for a single matrix or a pair of integer matrices
// the kappas are computed in a column-wise manner

// [[Rcpp::export]]

NumericMatrix kappaMtx(IntegerMatrix x, String method) {

  // result storage

  int n = x.ncol();

  IntegerMatrix index_pairs = intPairs(n);
  int n_pairs = index_pairs.nrow();

  NumericMatrix result(n_pairs, 4);

  colnames(result) =
    CharacterVector({"variable1", "variable2", "n", "kappa"});

  // computation of the kappas

  NumericVector pair_result(2);

  IntegerVector idx(2);

  for(int i = 0; i < n_pairs; ++i) {

    idx = index_pairs(i, _);

    pair_result = kappaCpp(x(_, idx[0]), x(_, idx[1]), method);

    result(i, 0) = 1.0 * idx[0] + 1.0;
    result(i, 1) = 1.0 * idx[1] + 1.0;
    result(i, 2) = pair_result[0];
    result(i, 3) = pair_result[1];

  }

  return result;

}

// [[Rcpp::export]]

NumericMatrix kappa2Mtx(IntegerMatrix x,
                        IntegerMatrix y,
                        String method) {

  /// the matrices need to have equal dimensions

  if((x.nrow() != y.nrow()) | (x.ncol() != y.ncol())) {

    stop("Matrices x and y must have equal dimensions");

  }

  /// computation

  int n_col = x.ncol();

  NumericMatrix res(n_col, 2);

  colnames(res) = CharacterVector({"n", "kappa"});

  for(int i = 0; i < n_col; ++i) {

    res(i, _) = kappaCpp(x(_, i), y(_, i), method);

  }

  return res;

}

// END
