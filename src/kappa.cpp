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

double kappaCpp(IntegerVector x,
                IntegerVector y,
                String  type) {

  /// contingency matrix: length checks and NA removal
  /// is done by the internal utility

  NumericMatrix ctg_mtx = irrTable(x, y);

  int n_col = ctg_mtx.ncol();

  if(n_col < 2) stop("At least two non-NA observations are required");

  /// matrix diagonal and margins

  NumericVector ctg_diag = mtxDiag(ctg_mtx);

  double ctg_total_freq = mtxTotalSum(ctg_mtx);

  NumericVector ctg_row_freq = mtxRowFreq(ctg_mtx);
  NumericVector ctg_col_freq = mtxColFreq(ctg_mtx);

  /// po: observed agreement (matrix diagonal by total frequency)
  /// pc: theoretical coincidence by chance
  /// used for unweighted kappa

  double po;
  double pc;

  if(type == "unweighted") {

    po = sum(ctg_diag)/ctg_total_freq;
    pc = crossProduct(ctg_col_freq, ctg_row_freq)(0, 0);

  } else {

    /// common variables for calculation of the weights

    NumericVector dim_seq = intSeq(1, n_col);

    NumericMatrix outer_delta = mtxAbs(outerDelta(dim_seq, dim_seq, false));

    outer_delta = mtxConstProd(outer_delta, 1.0/(n_col - 1));

    NumericMatrix weights;

    if(type == "equal") {

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

  return (po - pc)/(1 - pc);

}

// END

