/*** one-way ANOVA */

#include <Rcpp.h>
#include <Rmath.h>
#include "numericUtils.h"
#include "vectorUtils.h"
#include "transformUtils.h"
#include "oneAnova.h"

using namespace Rcpp;

// the base version of one-way ANOVA taking a list of numeric vectors

// [[Rcpp::export]]

NumericVector oneAnovaBase(List x) {

  // result container

  NumericVector result(11, NA_REAL);

  result[0] = 0.0;

  result.names() =
    CharacterVector::create("redFlag", "n", "k",
                            "ss_between", "ss_within", "ss_total",
                            "f", "df1", "df2",
                            "p_value", "etasq");

  // basic constants and containers

  int k = x.length(); // number of groups

  if(k < 2) {

    warning("At least two analysis groups are required");

    result[0] = 1.0;

    return result;

  }

  double grandMean = 0.0; // mean of all observations
  double N = 0.0; // will contain the complete observation number

  NumericVector groupSizes(k); //
  NumericVector groupMeans(k); // means of the groups

  // computation of group means and the grand mean

  for(int i = 0; i < k; ++i) {

    NumericVector group = x[i];

    group = na_omit(group);

    groupSizes[i] = group.length() * 1.0;

    if(groupSizes[i] < 2) {

      warning("Not enough observations to evaluate");

      result[0] = 1.0;

      return result;

    }

    N += groupSizes[i];
    grandMean += sum(group);

    groupMeans[i] = mean(group);

    x[i] = group;

  }

  grandMean = grandMean/N;

  // sum of squares

  double ss_between = 0.0;
  double ss_within = 0.0;
  double ss_total = 0.0;

  for(int i = 0; i < k; ++i) {

    NumericVector group = x[i];

    ss_between += groupSizes[i] * std::pow(groupMeans[i] - grandMean, 2.0);

    NumericVector withinDeltas = group - groupMeans[i];

    ss_within += sum(withinDeltas * withinDeltas);

    NumericVector totalDeltas = group - grandMean;

    ss_total += sum(totalDeltas * totalDeltas);

  }

  // degrees of freedom and the F statistic, p value, and // [[Rcpp::export]]
  // eta-square effect size

  double df1 = (k - 1) * 1.0; // nominator df
  double df2 = (N - k) * 1.0; // denominator df

  double F = (ss_between/df1)/(ss_within/df2);

  double p_value = R::pf(F, df1, df2, false, false);

  if(p_value > 1.0) p_value = 1.0;

  double etasq = ss_between/ss_total;

  // output

  result[1] = N;
  result[2] = k;

  result[3] = ss_between;
  result[4] = ss_within;
  result[5] = ss_total;

  result[6] = F;
  result[7] = df1;
  result[8] = df2;

  result[9] = p_value;
  result[10] = etasq;

  return result;

}

// a version for a numeric vector and a splitting factor

// [[Rcpp::export]]

NumericVector oneAnovaVec(NumericVector x,
                          IntegerVector f,
                          bool crash = false) {

  // processing of the vectors

  LogicalVector idx = !((is_na(x)) | (is_na(f)));

  x = x[idx];
  f = f[idx];

  List x_splits = Split(x, f);

  // testing

  NumericVector result = oneAnovaBase(x_splits);

  double redFlag = result[0];
  result = result[Range(1, 10)];

  if((redFlag == 1.0) & crash) stop("Computation failed");

  if(redFlag == 1.0) return result;

  result.names() =
    CharacterVector::create("n", "k",
                            "ss_between", "ss_within", "ss_total",
                            "f", "df1", "df2",
                            "p_value", "etasq");

  return result;

}

// and a version for numeric matrices

// [[Rcpp::export]]

NumericMatrix oneAnovaMtx(NumericMatrix x,
                          IntegerVector f,
                          bool crash = true) {

  int nVars = x.ncol();

  NumericMatrix result(nVars, 10);

  colnames(result) =
    CharacterVector::create("n", "k",
                            "ss_between", "ss_within", "ss_total",
                            "f", "df1", "df2",
                            "p_value", "etasq");

  NumericVector tstOutput(10);

  for(int i = 0; i < nVars; ++i) {

    tstOutput = oneAnovaVec(x(_, i), f, crash);

    result(i, _) = tstOutput;

  }

  if(checkNames(x)[1]) rownames(result) = colnames(x);

  return result;

}

// END
