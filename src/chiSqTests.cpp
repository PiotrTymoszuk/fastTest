/*** Tests with contingency tables */

#include <Rcpp.h>
#include <Rmath.h>
#include "contingencyUtils.h"
#include "chiSqTests.h"
#include "transformUtils.h"

using namespace Rcpp;

// checking the contingency table for possible flaws. The function returns
// true if any are found and prints a message.

bool ctgFlaws(NumericMatrix ctg) {

  bool res = false;
  int nRows = ctg.nrow();
  int nCols = ctg.ncol();

  // checking dimensions of the table

  if((nRows == 1) | (nCols == 1)) {

    res = true;

    Rcout << "One-dimensional contingency table" << "\n";

    return res;

  }

  // resolving the matrix content into a numeric vector

  NumericVector ctgVals(nRows * nCols);

  for(int i = 0; i < nRows; ++i) {

    for(int j = 0; j < ctg.ncol(); ++j) {

      ctgVals[i + j] = ctg(i, j);

    }

  }

  // checking if any values are negative or all of them are 0

  for(int i = 0; i < ctgVals.length(); ++i) {

    if(ctgVals[i] < 0) {

      res = true;

      Rcout << "There are negative values in the contingency matrix" << "\n";

      break;

    }

  }

  if(sum(ctgVals) < 1) {

    res = true;

    Rcout << "All values of the contingency matrix are zero" << "\n";

  }

  return res;

}

// Equivalent of R's chisq-test in Rcpp taking a contingency table
// as an argument

// [[Rcpp::export]]

NumericVector chiSqTestTbl(NumericMatrix ctg,
                           bool correct = true,
                           bool crash = true) {

  // checking the contingency matrix for flaws

  NumericVector res = NumericVector::create(Named("n") = NA_REAL,
                                            Named("chisq") = NA_REAL,
                                            Named("df") = NA_REAL,
                                            Named("p_value") = NA_REAL,
                                            Named("cramer_v") = NA_REAL);

  bool redFlag = ctgFlaws(ctg);

  if(redFlag & crash) stop("Critical contingency matrix errors");

  if(redFlag) {

    warning("Contingency matrix issues");

    return res;

  }

  // Compute row and column sums

  int nRows = ctg.nrow();
  int nCols = ctg.ncol();

  NumericVector rowSums(nRows);
  NumericVector colSums(nCols);

  for (int i = 0; i < nRows; ++i) {

    for (int j = 0; j < nCols; ++j) {

      rowSums[i] += ctg(i, j);
      colSums[j] += ctg(i, j);

    }

  }

  // Compute expected frequencies

  NumericMatrix expected(nRows, nCols);

  for (int i = 0; i < nRows; ++i) {

    for (int j = 0; j < nCols; ++j) {

      expected(i, j) = (rowSums[i] * colSums[j]) / sum(rowSums);

    }

  }

  // Compute chi-squared statistic

  double chiSquared = 0.0;
  double yates = 0.0;

  if(correct) yates = 0.5;

  for (int i = 0; i < nRows; ++i) {

    for (int j = 0; j < nCols; ++j) {

      double diff = ctg(i, j) - expected(i, j);

      if(yates < diff) diff = diff - yates;

      chiSquared += std::pow(diff, 2.0) / expected(i, j);

    }

  }

  // Compute degrees of freedom

  int df = (nRows - 1) * (nCols - 1);

  // Compute p-value

  double pValue = R::pchisq(chiSquared, df, false, false);

  // Compute Cramer's V effect size

  double denom = nRows - 1;

  if(nCols < nRows) denom = nCols - 1;
  if(denom == 0) denom = 1;

  double cramerV = std::pow(chiSquared/sum(rowSums), 0.5)/denom;

  res[0] = sum(rowSums);
  res[1] = chiSquared;
  res[2] = df;
  res[3] = pValue;
  res[4] = cramerV;

  return res;

}

// chi-square test for a numeric vector and a splitting integer variable

// [[Rcpp::export]]

NumericVector chiSqTestVec(NumericVector x,
                           IntegerVector f,
                           bool correct = true,
                           bool crash = true) {

  // contingency matrix, NAs are silently removed

  NumericMatrix ctg_mtx = xTable(x, f);

  return chiSqTestTbl(ctg_mtx, correct, crash);

}

// chi-square test for a numeric matrix and a splitting integer variable

// [[Rcpp::export]]

NumericMatrix chiSqTestMtx(NumericMatrix x,
                           IntegerVector f,
                           bool correct = true,
                           bool crash = true) {

  // testing

  int xColSize = x.ncol();

  NumericMatrix resMtx(xColSize, 5);

  for(int i = 0; i < xColSize; ++i) {

    NumericVector testInput = x(_, i);

    NumericVector testRes = chiSqTestVec(testInput, f, correct, crash);

    resMtx(i, _) = testRes;

  }

  // names and output

  colnames(resMtx) =
    CharacterVector::create("n", "chisq", "df", "p_value", "cramer_v");

  LogicalVector nameCheck = checkNames(x);

  if(nameCheck[1]) rownames(resMtx) = colnames(x);

  return resMtx;

}

// END
