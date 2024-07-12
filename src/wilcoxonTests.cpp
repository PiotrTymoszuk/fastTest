/*** Wilcoxon tests and effect sizes */

#include <Rcpp.h>
#include <Rmath.h>
#include "contingencyUtils.h"
#include "rankUtils.h"
#include "numericUtils.h"
#include "transformUtils.h"
#include "vectorUtils.h"
#include "wilcoxonTests.h"

using namespace Rcpp;

// biseral r correlation or effect size statistic for a numeric vector x,
// and a splitting integer vector f.
// This function computes the key stats for the Wilcoxon tests as well, i.e
// U1, U2, numbers of pairs in favor and against the alternative hypothesis,
// as well as numbers of observations or ranks not affected by ties.
//
// Additionally, the function performs quality control checks; the first element
// of the output numeric vector is a flag (0: no issues, 1: some issues present).
//
// The output vector positions:
// 0: issue flag,
// 1: numbers of complete observations in the first group,
// 2: number of observations in the second analysis group,
// 3: numbers of observation pairs without zero differences (NA for unpaired),
// 4: tie correction factor (`r sum(NTIES^3 -NTIES)/48`),
// 5: U1,
// 6: U2,
// 7: the total number of evaluated pairs,
// 8: fraction of unfavorable pairs (u),
// 9: fraction of favorable pairs (f),
// 10: Hodges-Lehman estimate of location or median delta.
// 11: biserial r statistic.
//
// the argument hl_estimate allows for computation of Hedges-Lehman estimates
// which may be however lengthy for large samples!!!

NumericVector biserialR(NumericVector x,
                        IntegerVector f,
                        bool paired = false,
                        bool hl_estimate = false) {

  // result container

  NumericVector result(12, NA_REAL);

  result.names() =
    CharacterVector::create("redFlag",
                            "n1", "n2",
                            "n_nonzero", "tie_correction",
                            "U1", "U2", "pair_number",
                            "h0_fraction", "h1_fraction",
                            "estimate", "biserial_r");

  result[0] = 0;

  // general entry control

  IntegerVector f_cats = unique(f);
  f_cats = na_omit(f_cats);

  if(f_cats.length() < 2) {

    warning("Only one strata present");

    result[0] = 1;

    return result;

  }

  if(f_cats.length() > 2) {

    warning("More than two strata present");

    result[0] = 1;

    return result;

  }

  // Mann-Whitney statistic containers, numbers of untied pairs,
  // numbers of complete observations

  double S;

  double U1;
  double U2;

  int n; // untied pair number: paired test
  int n1; // numbers of complete observations in the first group
  int n2; // numbers of complete observations in the second group

  double delta_est; // difference of medians
  NumericVector ranks; //stores the ranks or signed ranks

  if(paired) {

    // checking for the proper strata sizes and selecting complete cases

    List x_splits = Split(x, f);

    NumericVector x1 = x_splits[0];
    NumericVector x2 = x_splits[1];

    List procRes = procVectors(x1, x2, true);

    x1 = procRes[0];
    x2 = procRes[1];

    if(procRes[2]) result[0] = 1;

    // signed ranks for differences between the strata

    NumericVector vec_diff = delta(x1, x2);

    NumericVector vec_ranks = signedRank(vec_diff, "average", false);

    n1 = x1.length();
    n2 = x2.length();

    n = vec_ranks.length();

    if(n  < 2) {

      warning("Not enough untied observations to compute the test.");

      result[0] = 1;

      return result;

    }

    S = (n * (n + 1)/2);

    // finding and summing ranks for pairs with positive and negative
    // differences

    LogicalVector index1 = vec_ranks > 0;
    LogicalVector index2 = vec_ranks < 0;

    NumericVector ranks1 = vec_ranks[index1];
    NumericVector ranks2 = vec_ranks[index2];

    // the swap below: to keep the level order (level 2 against level 1)

    U2 = sum(ranks1);
    U1 = sum(ranks2) * -1;

    result[3] = n;

    if(hl_estimate) {

      delta_est = locationPair(vec_diff); // single sample Hodges-Lehman

    } else {

      delta_est = Median(vec_diff, true);

    }

    ranks = Abs(vec_ranks, true); // to keep it consistent with the R's wilcox.test

  } else {

    // complete cases

    List compCases = completeCases(x, f);

    x = compCases[0];
    f = compCases[1];

    // setting common ranks for the x vector and splitting by
    // the vector f

    NumericVector x_ranks = Rank(x, "average", false);

    List rank_splits = Split(x_ranks, f);

    NumericVector ranks1 = rank_splits[0];
    NumericVector ranks2 = rank_splits[1];

    ranks = x_ranks;

    n1 = ranks1.length();
    n2 = ranks2.length();

    if((n1 < 2) | (n2 < 2)) {

      warning("Not enough complete observations to compute the test");

      result[0] = 1;

      return result;

    }

    S = n1 * n2;

    U1 = sum(ranks1) - n1 * (n1 + 1)/2;
    U2 = sum(ranks2) - n2 * (n2 + 1)/2;

    // Hodges-Lehman estimates

    List x_splits = Split(x, f);

    if(hl_estimate) {

      delta_est = locationStd(x_splits[0], x_splits[1]);

    } else {

      double med1 = Median(x_splits[0], true);
      double med2 = Median(x_splits[1], true);

      delta_est = med2 - med1;

    }

  }

  // frequency of pairs in favor of the alternative hypothesis (f_pairs)
  // and in favor of the NULL hypothesis (u_pairs)

  double u_pairs = U1/S;
  double f_pairs = U2/S;

  double rbs = f_pairs - u_pairs;

  // calculation of numbers of ties and tie corrections

  IntegerVector rank_tbl = Table(ranks);

  IntegerVector tieCorrection = rank_tbl * rank_tbl * rank_tbl - rank_tbl;

  result[4] = sum(tieCorrection);

  // the output vector

  result[1] = n1;
  result[2] = n2;

  result[5] = U1;
  result[6] = U2;
  result[7] = S;
  result[8] = u_pairs;
  result[9] = f_pairs;
  result[10] = delta_est;
  result[11] = rbs;

  return result;

}

// Wilcoxon test return the following stats:
// 0: number of cases in the first group,
// 1: number of cases in the second group,
// 2: U statistic,
// 3: p value
// 4: total number of evaluated pairs, i.e. S,
// 5: fraction of pairs against the alternative hypothesis,
// 6: number of fractions in favor of the alternative hypothesis,
// 7: estimate: HL estimate or differerence of medians
// 8: biserial r as an effect size statistic
//
// core variant of the Wilcoxon test working with vectors

// [[Rcpp::export]]

NumericVector testWilcoxonVec(NumericVector x,
                              IntegerVector f,
                              String type = "standard",
                              String alternative = "two.sided",
                              bool exact = false,
                              bool correct = true,
                              bool hl_estimate = false,
                              bool crash = true) {

  // computation of the N numbers, U stats, pairs and biseral r, as well as
  // a check for possible flaws

  bool paired = false;

  if(type == "paired") paired = true;

  NumericVector result(9, NA_REAL);

  result.names() =
    CharacterVector::create("n1", "n2", "u", "p_value",
                            "pair_number", "h0_fraction", "h1_fraction",
                            "estimate", "biserial_r");

  NumericVector stats = biserialR(x, f, paired, hl_estimate);

  result[0] = stats[1];
  result[1] = stats[2];

  result[4] = stats[7];
  result[5] = stats[8];
  result[6] = stats[9];

  result[7] = stats[10];
  result[8] = stats[11];

  if(crash & (stats[0] == 1)) stop("Critical input vector issues");

  if(stats[0] == 1) {

    warning("Input vector issues");

    return result;

  }

  // test stat container

  double U; // Mann-Whitney U statistic
  double n; // number of samples
  double tie_correction = stats[4]; // number of ties
  double z = 1.0; // normal distribution z statistic: normal approximation
  double sigma = 1.0; // sigma of the normal distribution

  double p_value = 1.0;

  // testing

  if(paired) {

    // wilcoxon's signed rank test

    U = stats[6];
    n = stats[3];

    double n_zeros = stats[2] - stats[3];

    if(exact & (tie_correction == 0) & (n_zeros == 0)) {

      // exact p value

      if(alternative == "two.sided") {

        if (U > (n * (n + 1)/4)) {

          p_value = R::psignrank(U - 1, n, false, false);

        } else {

          p_value = R::psignrank(U, n, true, false);

        }

        p_value = 2 * p_value;

      } else if(alternative == "greater") {

        p_value = R::psignrank(U - 1, n, false, false);

      } else {

        p_value = R::psignrank(U, n, true, false);

      }

    } else {

      // normal distribution approximation

      z = U - n * (n + 1)/4;

      sigma = std::sqrt(n * (n + 1) * (2 * n + 1)/24 - tie_correction/48);

      double correction = 0.0; // continuity correction, specified later

      if(alternative == "less") {

        if(correct) correction = -0.5;

        z = (z - correction)/sigma;

        p_value = R::pnorm(z, 0.0, 1.0, true, false);

      } else if(alternative == "greater") {

        if(correct) correction = 0.5;

        z = (z - correction)/sigma;

        p_value = R::pnorm(z, 0.0, 1.0, false, false);

      } else {

        if(correct) {

          if(z == 0) {

            correction = 0.0;

          } else {

            correction = 0.5 * z/std::abs(z);

          }

        }

        z = (z - correction)/sigma;

        p_value = R::pnorm(z, 0.0, 1.0, true, false);
        double p2 = R::pnorm(z, 0.0, 1.0, false, false);

        if(p2 < p_value) p_value = p2;

        p_value = 2 * p_value;

      }

    }

  } else {

    // Wilcoxon's rank test also known as Mann-Whitney U test

    U = stats[6];

    if(stats[5] < U) U = stats[5];

    double n_x = stats[1];
    double n_y = stats[2];

    // defining the W statistic and using it in the code below: R legacy
    // in out case this is U2, i.e. U statistic for the second analysis group

    double W = stats[6];

    if(exact & (tie_correction == 0)) {

      // exact variant, p value derived from the Wilcoxon's distribution
      // only in case there are no ties

      if(alternative == "two.sided") {

        if(W > (n_x * n_y/2)) {

          p_value = R::pwilcox(W - 1, n_y, n_x, false, false);

        } else {

          p_value = R::pwilcox(W, n_y, n_x, true, false);

        }

        p_value = p_value * 2;

      } else if(alternative == "greater") {

        p_value = R::pwilcox(W - 1, n_y, n_x, false, false);

      } else {

        p_value = R::pwilcox(W, n_y, n_x, true, false);

      }

    } else {

      z = W - n_x * n_y/2;

      sigma =
        std::sqrt((n_x * n_y/12) * ((n_x + n_y + 1) - tie_correction/((n_x + n_y) * (n_x + n_y - 1))));

      double correction = 0.0; // continuity correction, specified later

      if(alternative == "less") {

        if(correct) correction = -0.5;

        z = (z - correction)/sigma;

        p_value = R::pnorm(z, 0, 1, true, false);

      } else if(alternative == "greater") {

        if(correct) correction = 0.5;

        z = (z - correction)/sigma;

        p_value = R::pnorm(z, 0, 1, false, false);

      } else {

        if(correct) {

          if(z == 0) {

            correction = 0;

          } else {

            correction = z/std::abs(z) * 0.5;

          }

        }

        z = (z - correction)/sigma;

        p_value = R::pnorm(z, 0, 1, true, false);

        double p2 = R::pnorm(z, 0, 1, false, false);

        if(p2 < p_value) p_value = p2;

        p_value = p_value * 2;

      }

    }

  }

  // output

  result[2] = U;

  if(p_value > 1) p_value = 1;

  result[3] = p_value;

  return result;

}

// variant of the function taking two numeric vectors

// [[Rcpp::export]]

NumericVector testWilcoxonStd(NumericVector x,
                              NumericVector y,
                              String type = "standard",
                              String alternative = "two.sided",
                              bool exact = false,
                              bool correct = true,
                              bool hl_estimate = false,
                              bool crash = true) {

  // preparing the input for `testWilcoxonVec()`

  NumericVector newX = Concat(x, y);
  NumericVector numF = Concat(rep(1.0, x.length()), rep(2.0, y.length()));

  IntegerVector f = as<IntegerVector>(numF);

  return testWilcoxonVec(newX, f, type, alternative, exact, correct, hl_estimate, crash);

}

// variant of the function working

// [[Rcpp::export]]

NumericMatrix testWilcoxonMtx(NumericMatrix x,
                              IntegerVector f,
                              String type = "standard",
                              String alternative = "two.sided",
                              bool exact = false,
                              bool correct = true,
                              bool hl_estimate = false,
                              bool crash = true) {

  // entry control and containers

  int xColSize = x.ncol();

  NumericMatrix resMtx(xColSize, 9);

  for(int i = 0; i < xColSize; ++i) {

    NumericVector testInput = x(_, i);

    NumericVector testRes =
      testWilcoxonVec(testInput, f, type, alternative, exact, correct, hl_estimate, crash);

    resMtx(i, _) = testRes;

  }

  // naming and output

  colnames(resMtx) =
    CharacterVector::create("n1", "n2", "u", "p_value",
                            "pair_number", "h0_fraction", "h1_fraction",
                            "estimate", "biserial_r");

  LogicalVector nameCheck = checkNames(x);

  if(nameCheck[1]) rownames(resMtx) = colnames(x);

  return resMtx;

}

// END
