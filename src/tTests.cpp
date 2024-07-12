/*** Cohen's d effect size and T tests */

#include <Rcpp.h>
#include <Rmath.h>
#include "numericUtils.h"
#include "contingencyUtils.h"
#include "vectorUtils.h"
#include "transformUtils.h"
#include "tTests.h"

using namespace Rcpp;

// pooled variance and Cohen's effect size

NumericVector cohenD(NumericVector x, NumericVector y, String type = "standard") {

  double (*fun)(NumericVector, NumericVector);

  // computing the pooled variance

  if(type == "standard") {

    fun = poolVarStandard;

  } else if(type == "welch") {

    fun = poolVarWelch;

  } else if(type == "paired") {

    fun = poolVarPaired;

  }

  double poolVar = fun(x, y);

  if(poolVar == 0) {

    warning("Pooled variance is zero: identical input vectors?");

    return NumericVector::create(0, NA_REAL);

  }

  double mean_diff;

  if(type == "paired") {

    NumericVector vec_diff = delta(x, y);

    mean_diff = mean(vec_diff);

  } else {

    mean_diff = mean(y) - mean(x);

  }

  return NumericVector::create(poolVar, mean_diff/std::sqrt(poolVar));

}

// T tests return the following:
// 0: number of cases in the first group, 1: number of cases in the second group,
// 2: T statistic, 3: number of degrees of freedom, 4: p-value,
// 5: mean in the first group, 6: mean in the second group, 7: difference of
// means, 8: lower 95% CI bound, 9: upper 95% CI bound,
// 10: pooled variance, 11: Cohen's d

// T test operating with two numeric vectors

// [[Rcpp::export]]

NumericVector tTestStd(NumericVector x,
                       NumericVector y,
                       String type = "standard",
                       String alternative = "two.sided",
                       double conf_level = 0.95,
                       bool crash = true) {

  // check for the vectors

  NumericVector testRes(13, NA_REAL);

  testRes.names() =
    CharacterVector::create("n1", "n2", "t", "df", "p_value",
                            "mean1", "mean2", "estimate", "sem",
                            "lower_ci", "upper_ci",
                            "pooled_var", "cohen_d");

  bool paired = false;

  if(type == "paired") paired = true;

  List vectorProc = procVectors(x, y, paired);

  bool redFlags = vectorProc[2];

  if(crash & redFlags) stop("Critical issues of the input vectors");

  if(redFlags) {

    warning("Issues of the input vectors");

    return testRes;

  }

  x = vectorProc[0];
  y = vectorProc[1];

  // vector sizes, means, and variances

  int xLen = x.length();
  int yLen = y.length();

  double mean_x = mean(x);
  double mean_y = mean(y);

  double mean_diff = mean_y - mean_x;

  double var_x = Var(x, true);
  double var_y = Var(y, true);

  double df = 1.0;
  double pooled_var = 1.0;
  double pooled_sem = 1.0;
  double t_stat = 1.0;

  // pooled variance, degrees of freedom, and T statistic

  if(type == "standard") {

    // Student's T test variant

    df = xLen + yLen - 2;

    pooled_var = ((xLen - 1) * var_x + (yLen - 1) * var_y) / df;
    pooled_sem = std::sqrt(pooled_var * (1.0/xLen + 1.0/yLen));

  } else if(type == "welch") {

    // Welch's T test variant

    double sem_x = std::sqrt(var_x/xLen);
    double sem_y = std::sqrt(var_y/yLen);

    pooled_sem = sqrt(std::pow(sem_x, 2.0) + std::pow(sem_y, 2.0));

    double df_nom = std::pow(pooled_sem, 4.0);
    double df_denom = std::pow(sem_x, 4.0)/(xLen - 1) + std::pow(sem_y, 4.0)/(yLen - 1);

    df = df_nom/df_denom;

  } else {

    // paired T test

    NumericVector vec_diff = delta(x, y);
    mean_x = NA_REAL;
    mean_y = NA_REAL;

    pooled_sem = std::sqrt(Var(vec_diff, true)/xLen);
    mean_diff = mean(vec_diff);
    df = xLen - 1;

  }

  t_stat = mean_diff/pooled_sem;

  // p value, degrees of freedom, and confidence intervals

  double p_value = 1.0;
  NumericVector ci(2, NA_REAL);

  if(alternative == "less") {

    p_value = R::pt(-t_stat, df, true, false);

    ci[0] = R_NegInf;
    ci[1] = t_stat + R::qt(conf_level, df, true, false);

  } else if(alternative == "greater") {

    p_value = R::pt(t_stat, df, true, false);

    ci[0] = t_stat - R::qt(conf_level, df, true, false);
    ci[1] = R_PosInf;

  } else {

    p_value = 2.0 * R::pt(std::abs(t_stat), df, false, false);

    double alpha = 1 - conf_level;

    ci[1] = R::qt(1 - alpha/2, df, true, false);
    ci[0] = -ci[1];

    ci[0] = t_stat + ci[0];
    ci[1] = t_stat +  ci[1];

  }

  ci[0] = ci[0] * pooled_sem;
  ci[1] = ci[1] * pooled_sem;

  // effect size

  NumericVector effSize = cohenD(x, y, type);

  testRes[0] = xLen;
  testRes[1] = yLen;
  testRes[2] = t_stat;
  testRes[3] = df;
  testRes[4] = p_value;
  testRes[5] = mean_x;
  testRes[6] = mean_y;
  testRes[7] = mean_diff;
  testRes[8] = pooled_sem;
  testRes[9] = ci[0];
  testRes[10] = ci[1];
  testRes[11] = effSize[0];
  testRes[12] = effSize[1];

  return testRes;

}

// T test operating with a numeric vector and a splitting factor

//[[Rcpp::export]]

NumericVector tTestVec(NumericVector x,
                       IntegerVector f,
                       String type = "standard",
                       String alternative = "two.sided",
                       double conf_level = 0.95,
                       bool crash = true) {

  List x_splits = Split(x, f);

  if(x_splits.length() < 2) stop("Only one strata available");

  return tTestStd(x_splits[0], x_splits[1], type, alternative, conf_level, crash);

}

// T test operating with a numeric matrix

//[[Rcpp::export]]

NumericMatrix tTestMtx(NumericMatrix x,
                       IntegerVector f,
                       String type = "standard",
                       String alternative = "two.sided",
                       double conf_level = 0.95,
                       bool crash = true) {

  // testing

  int xColSize = x.ncol();

  NumericMatrix resMtx(xColSize, 13);

  for(int i = 0; i < xColSize; ++i) {

    NumericVector testInput = x(_, i);

    NumericVector testRes =
      tTestVec(testInput, f, type, alternative, conf_level, crash);

    resMtx(i, _) = testRes;

  }

  // naming and output

  colnames(resMtx) =
    CharacterVector::create("n1", "n2", "t", "df", "p_value",
                            "mean1", "mean2", "estimate", "sem",
                            "lower_ci", "upper_ci",
                            "pooled_var", "cohen_d");

  LogicalVector nameCheck = checkNames(x);

  if(nameCheck[1]) rownames(resMtx) = colnames(x);

  return resMtx;

}

// END
