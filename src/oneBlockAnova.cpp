/*** one-way ANOVA for block designs with a single blocking factor */

#include <Rcpp.h>
#include <Rmath.h>
#include "numericUtils.h"
#include "vectorUtils.h"
#include "transformUtils.h"
#include "oneBlockAnova.h"

using namespace Rcpp;

// a base function for a numeric vector x, the treatment split factor f,
// and the block split factor b

// [[Rcpp::export]]

NumericVector oneBlockAnovaVec(NumericVector x,
                               IntegerVector f,
                               IntegerVector b) {

  // the result container

  NumericVector result(18, NA_REAL);

  result[0] = 0.0;

  result.names() =
    CharacterVector::create("redFlag", "n", "k", "b",
                            "ss_treatement", "ss_block", "ss_within", "ss_total",
                            "f_treatment", "df1_treatment", "df2_treatment",
                            "p_treatment", "etasq_treatment",
                            "f_block", "df1_block", "df2_block",
                            "p_block", "etasq_block");

  // cNA removal

  LogicalVector idx = !(is_na(x) | is_na(f) | is_na(b));

  x = x[idx];
  f = f[idx];
  b = b[idx];

  // basic containers and input control

  int k = unique(f).length(); // number of treatment groups

  if(k < 2) {

    warning("At least two analysis groups are required");

    result[0] = 1.0;

    return result;

  }

  int B = unique(b).length();

  if(B < 2) {

    warning("At least two observation blocks are required");

    result[0] = 1.0;

    return result;

  }

  bool fEqual = checkEqualSplits(as<NumericVector>(f));
  bool bEqual = checkEqualSplits(as<NumericVector>(b));

  if(!fEqual | !bEqual) {

    warning("Level-deficient design");

    result[0] = 1.0;

    return result;

  }

  // storage of means and sizes

  double grandMean = mean(x); // the grand mean
  double N = x.length(); // the total number of complete observations

  NumericVector groupMeans(k); // will store means of the treatment groups
  NumericVector groupSizes(k); // will store sizes of the treatment groups
  double ss_between = 0.0; // will store sum of squares for the treatment groups

  NumericVector blockMeans(B); // will store means of the blocks
  NumericVector blockSizes(B); // will store sizes of the blocks
  double ss_block = 0.0; // will store sum of squares for the blocks

  // computation of the group sizes and means

  List group_splits = Split(x, f);

  for(int i = 0; i < k; ++i) {

    NumericVector group = group_splits[i];

    if(group.length() < 2) {

      warning("Not enough observations in one of the treatment groups");

      result[0] = 1.0;

      return result;

    }

    groupSizes[i] = group.length() * 1.0;

    double groupAvg = mean(group);

    groupMeans[i] = groupAvg;

    ss_between += std::pow(groupAvg - grandMean, 2.0);

  }

  ss_between = B * ss_between;

  // calculation of the block sizes and means

  List block_splits = Split(x, b);

  for(int i = 0; i < B; ++i) {

    NumericVector block = block_splits[i];

    if(block.length() < 2) {

      warning("Not enough observations in one of the blocks");

      result[0] = 1.0;

      return result;

    }

    blockSizes[i] = block.length();

    double blockAvg = mean(block);

    blockMeans[i] = blockAvg;

    ss_block += std::pow(blockAvg - grandMean, 2.0);

  }

  ss_block = k * ss_block;

  // the total sum of squares,
  // within sum of squares will be calculated as a difference
  // of total sum of squares - between sum of squares - block sum of squares

  double ss_total = 0.0;

  for(int i = 0; i < N; ++i) ss_total += std::pow(x[i] - grandMean, 2.0);

  double ss_within = ss_total - ss_between - ss_block;

  // degrees of freedom and F stats, effect sizes

  double df1_treatment = k - 1; // nominator df
  double df2_treatment = (k - 1) * (B - 1); // denominator df for the within effects

  double F_treatment = (ss_between/df1_treatment)/(ss_within/df2_treatment);

  double etasq_treatment = ss_between/ss_total;


  double df1_block = B - 1;
  double df2_block = (k - 1) * (B - 1); // denominator df for the within effects

  double F_block = (ss_block/df1_block)/(ss_within/df2_treatment);

  double etasq_block = ss_block/ss_total;

  // p values from F distribution for the treatment effect

  double p2;

  double p_treatment = R::pf(F_treatment, df1_treatment, df2_treatment, true, false);

  p2 = R::pf(F_treatment, df1_treatment, df2_treatment, false, false);

  if(p2 < p_treatment) p_treatment = p2;
  if(p_treatment > 1.0) p_treatment = 1.0;

  // p values for the block effect

  double p_block = R::pf(F_block, df1_block, df2_block, true, false);

  p2 = R::pf(F_block, df1_block, df2_block, false, false);

  if(p2 < p_block) p_block = p2;
  if(p_block > 1.0) p_block = 1.0;

  // the output

  result[1] = N;
  result[2] = k;
  result[3] = B;

  result[4] = ss_between;
  result[5] = ss_block;
  result[6] = ss_within;
  result[7] = ss_total;

  result[8] = F_treatment;
  result[9] = df1_treatment;
  result[10] = df2_treatment;
  result[11] = p_treatment;
  result[12] = etasq_treatment;

  result[13] = F_block;
  result[14] = df1_block;
  result[15] = df2_block;
  result[16] = p_block;
  result[17] = etasq_block;

  return result;

}

// a matrix version

// [[Rcpp::export]]

NumericMatrix oneBlockAnovaMtx(NumericMatrix x,
                               IntegerVector f,
                               IntegerVector b,
                               bool crash = true) {

  // result container

  int nVars = x.ncol();

  NumericMatrix result(nVars, 17);

  colnames(result) =
    CharacterVector::create("n", "k", "b",
                            "ss_treatement", "ss_block", "ss_within", "ss_total",
                            "f_treatment", "df1_treatment", "df2_treatment",
                            "p_treatment", "etasq_treatment",
                            "f_block", "df1_block", "df2_block",
                            "p_block", "etasq_block");

  NumericVector tstResult(18);

  for(int i = 0; i < nVars; ++i) {

    tstResult = oneBlockAnovaVec(x(_, i), f, b);

    if((tstResult[0] == 1.0) & crash) stop("Calculation failed.");

    result(i, _) = tstResult[Range(1, 17)];

  }

  if(checkNames(x)[1]) rownames(result) = colnames(x);

  return result;

}

// END
