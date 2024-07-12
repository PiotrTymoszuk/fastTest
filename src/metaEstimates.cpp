/*** Meta-estimates caclulated with inverse variance method */

#include <Rcpp.h>
#include <Rmath.h>
#include "contingencyUtils.h"
#include "rankUtils.h"
#include "numericUtils.h"
#include "transformUtils.h"
#include "vectorUtils.h"
#include "metaEstimates.h"

using namespace Rcpp;

// base, fixed effect inverse variance meta-estimate
// https://meta-analysis.com/download/M-a_f_e_v_r_e_sv.pdf

NumericVector MetaFun(NumericVector y,
                      NumericVector e,
                      double tausq = 0.0,
                      String alternative = "two.sided",
                      double conf_level = 0.95) {

  // containers and quality control

  NumericVector result(6, NA_REAL);

  result.names() =
    CharacterVector::create("n",
                            "estimate", "sem",
                            "lower_ci", "upper_ci",
                            "z");

  int n = y.length();

  // common inverse variance and y-hat estimates

  NumericVector invVar(n, 1.0);
  NumericVector y_hat(n);

  for(int i = 0; i < n; ++i) {

    y_hat[i] = y[i]/(std::pow(e[i], 2.0) + tausq);

    invVar[i] = 1/(std::pow(e[i], 2.0) + tausq);

  }

  double denom = sum(invVar);

  double estimate = sum(y_hat)/denom;

  double sem = std::sqrt(1/denom);

  // z statistic, confidence intervals, and p value

  double z = estimate/sem;
  double alpha = 1 - conf_level;

  double lower_ci;
  double upper_ci;

  if(alternative == "two.sided") {

    alpha = alpha/2;

    lower_ci = estimate + R::qnorm(1 - alpha, 0, 1, false, false) * sem;
    upper_ci = estimate + R::qnorm(1 - alpha, 0, 1, true, false) * sem;

  } else if(alternative == "less") {

    lower_ci = R_NegInf;
    upper_ci = estimate + R::qnorm(conf_level, 0, 1, true, false) * sem;

  } else {

    lower_ci = estimate + R::qnorm(conf_level, 0, 1, false, false) * sem;
    upper_ci = R_PosInf;

  }

  // the output vector

  result[0] = n;
  result[1] = estimate;
  result[2] = sem;
  result[3] = lower_ci;
  result[4] = upper_ci;
  result[5] = z;

  return result;

}

// random effect variance estimation
// tau-square is estimated by the DerSimonian-Laird method

NumericVector estimateVar(NumericVector y, NumericVector e) {

  // calculation of the total variance Q, df, and tau-square: between-study
  // variance

  int n = y.length();

  NumericVector fixedWt(n);
  NumericVector q1(n);
  NumericVector q2(n);

  for(int i = 0; i < n; ++i) {

    fixedWt[i] = 1/std::pow(e[i], 2.0);

    q1[i] = fixedWt[i] * std::pow(y[i], 2.0);

    q2[i] = fixedWt[i] * y[i];

  }

  double Q = sum(q1) - std::pow(sum(q2), 2.0)/sum(fixedWt);

  double df = n - 1;

  double tausq = 0.0;

  double C = 1.0;

  if(Q > df) {

    C = sum(fixedWt) - sum(product(fixedWt, fixedWt))/sum(fixedWt);

    tausq = (Q - df)/C;

  }

  // output

  NumericVector result =
    NumericVector::create(Q, C, df, tausq);

  result.names() =
    CharacterVector::create("Q", "C", "df", "tausq");

  return result;

}

// the R-exposed Rcpp function for calculation of meta-estimated with the
// fixed method and the DerSimonian-Laird random-effect algorithm

//[[Rcpp::export]]

NumericVector metaVec(NumericVector y,
                      NumericVector e,
                      String type = "fixed",
                      String alternative = "two.sided",
                      double conf_level = 0.95,
                      bool crash = true) {

  // result container

  NumericVector result(10, NA_REAL);

  result.names() =
    CharacterVector::create("n", "q", "p_q_value", "tausq",
                            "estimate", "sem", "lower_ci", "upper_ci",
                            "z", "p_value");

  // entry control for the input vectors, NA removal

  List procVecs = procVectors(y, e, true);

  bool redFlags = procVecs[2];

  y = procVecs[0];
  e = procVecs[1];

  int n = y.length();

  for(int i = 0; i < n; ++i) {

    if(e[i] == 0.0) {

      redFlags = true;

      warning("Zeros detected in the error vector");

      break;

    }

  }

  if(redFlags & crash) stop("Critical issues of the input vectors");

  if(redFlags) {

    warning("Issues of the input vectors");

    return result;

  }

  // estimation of the between-study variance, i.e. Q and tau,
  // chi-square test for significance of the between-study effects

  NumericVector varStats = estimateVar(y, e);

  double p_q_value = R::pchisq(varStats[0], n - 1, false, false);

  // meta-estimates, fixed or random, as requested by the user

  NumericVector metaStats;

  if(type == "fixed") {

    metaStats = MetaFun(y, e, 0.0, alternative, conf_level);

  } else {

    metaStats = MetaFun(y, e, varStats[3], alternative, conf_level);

  }

  // p values of the pooled effects obtained from Z distribution

  double z = metaStats[5];
  double p_value;

  if(alternative == "less") {

    p_value = R::pnorm(z, 0.0, 1.0, true, false);

  } else if(alternative == "greater") {

    p_value = R::pnorm(z, 0.0, 1.0, false, false);

  } else {

    p_value = R::pnorm(z, 0.0, 1.0, true, false);
    double p2 = R::pnorm(z, 0.0, 1.0, false, false);

    if(p2 < p_value) p_value = p2;

    p_value = 2 * p_value;

  }

  if(p_value > 1) p_value = 1;

  // output

  result[0] = n;
  result[1] = varStats[0];
  result[2] = p_q_value;
  result[3] = varStats[3];

  result[4] = metaStats[1];
  result[5] = metaStats[2];
  result[6] = metaStats[3];
  result[7] = metaStats[4];

  result[8] = z;
  result[9] = p_value;

  return result;

}

// functions for calculation of meta-estimates with matrices and list as
// an input

// [[Rcpp::export]]

NumericMatrix metaMtx(NumericMatrix y,
                      NumericMatrix e,
                      String type = "fixed",
                      String alternative = "two.sided",
                      double conf_level = 0.95,
                      bool crash = true) {

  // entry control

  int n = y.ncol();

  if(n != e.ncol()) stop("Input matrices of non-equal dimensions");

  if(y.nrow() != e.nrow()) stop("Input matrices of non-equal dimensions");

  // result container

  NumericMatrix result(n, 10);

  colnames(result) =
    CharacterVector::create("n", "q", "p_q_value", "tausq",
                            "estimate", "sem", "lower_ci", "upper_ci",
                            "z", "p_value");

  if(checkNames(y)[1]) {

    rownames(result) = colnames(y);

  } else if(checkNames(y)[1]) {

    rownames(result) = colnames(y);

  }

  // calculation of the meta-estimates

  NumericVector inputY(y.nrow());
  NumericVector inputE(y.nrow());

  NumericVector metaStats(10, NA_REAL);

  for(int i = 0; i < n; ++i) {

    inputY = y(_, i);
    inputE = e(_, i);

    metaStats = metaVec(inputY, inputE, type, alternative, conf_level, crash);

    result(i, _) = metaStats;

  }

  return result;


}


// [[Rcpp::export]]

NumericMatrix metaList(List y,
                       List e,
                       String type = "fixed",
                       String alternative = "two.sided",
                       double conf_level = 0.95,
                       bool crash = true) {

  // input control

  int n = y.length();

  if(n != e.length()) stop("Incompatible lengths of input vectors");

  // result container

  NumericMatrix result(n, 10);

  colnames(result) =
    CharacterVector::create("n", "q", "p_q_value", "tausq",
                            "estimate", "sem", "lower_ci", "upper_ci",
                            "z", "p_value");

  if(checkListNames(y)) {

    CharacterVector new_names = y.names();

    rownames(result) = new_names;

  } else if(checkListNames(y)) {

    CharacterVector new_names = e.names();

    rownames(result) = new_names;

  }

  // calculation of the meta-estimates

  NumericVector metaStats(10, NA_REAL);

  for(int i = 0; i < n; ++i) {

    metaStats = metaVec(y[i], e[i], type, alternative, conf_level, crash);

    result(i, _) = metaStats;

  }

  return result;

}

// END
