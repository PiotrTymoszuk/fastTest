/*** Levene test-calculating functions */

#include <Rcpp.h>
#include <Rmath.h>
#include "numericUtils.h"
#include "vectorUtils.h"
#include "transformUtils.h"
#include "leveneTests.h"

using namespace Rcpp;

// Levene's test following the example from Wikipedia
// https://en.wikipedia.org/wiki/Levene%27s_test
// implementing the Brown-Forsythe variant, i.e. with group medians instead of
// group means used for calculation of the centrality statistic of the group
// and the scores Z.

// base version for a list

// [[Rcpp::export]]

NumericVector leveneBase(List x, String type = "standard") {

  // result container

  NumericVector result(7);

  result.names() =
    CharacterVector::create("redFlag", "n", "k", "f", "df1", "df2", "p_value");

  result[0] = 0.0;

  // entry control and base containers, calculation of the Z-scores

  int k = x.length(); // number of groups

  if(k < 2) {

    warning("At least two groups are required");

    result[0] = 1.0;

    return(result);

  }

  int N = 0; // total number of cases
  IntegerVector groupN(k); // vector of numbers of cases in each group

  List zList(k); // a list to store the Z-scores

  double grandMean = 0.0; // will store the grand mean of Z scores
  NumericVector zMeans(k); // will store means of Z-scores in the groups

  for(int i = 0; i < k; ++i) {

    NumericVector group = x[i];
    group = na_omit(group);

    groupN[i] = group.length();

    if(groupN[i] < 2) {

      warning("Not enough observations to evaluate: n < 2");

      result[0] = 1.0;

      return result;

    }

    N += group.length();

    double groupCenter = 0.0;

    if(type == "standard") {

      groupCenter = mean(group);

    } else {

      groupCenter = Median(group, false);

    }

    NumericVector zScores = Abs(group - groupCenter, false);
    zList[i] = zScores;

    grandMean += sum(zScores);
    zMeans[i] = mean(zScores);

  }

  grandMean = grandMean/N;

  // calculation of the test statistic F

  double df1 = 1.0 * (k - 1);
  double df2 = 1.0 * (N - k);

  NumericVector ssNom(k); // will store squares of between-group differences
  NumericVector ssDenom(k); // will store squares of within-group differences

  for(int i = 0; i < k; ++i) {

    ssNom[i] = (groupN[i] * 1.0) * std::pow((zMeans[i] - grandMean), 2.0);

    NumericVector zScores = zList[i];
    NumericVector zDeltas = zScores - zMeans[i];

    ssDenom[i] = sum(zDeltas * zDeltas);

  }

  double F = df2/df1 * sum(ssNom)/sum(ssDenom);

  // calculation of p values and output

  double p_value = R::pf(F, df1, df2, false, false);

  if(p_value > 1.0) p_value = 1.0;

  result[1] = N;
  result[2] = k;
  result[3] = F;
  result[4] = df1;
  result[5] = df2;
  result[6] = p_value;

  return result;

}

// version for a numeric vector and a splitting factor

// [[Rcpp::export]]

NumericVector leveneVec(NumericVector x,
                        IntegerVector f,
                        String type = "standard",
                        bool crash = true) {

  // processing of the vectors

  LogicalVector idx = !((is_na(x)) | (is_na(f)));

  x = x[idx];
  f = f[idx];

  List x_splits = Split(x, f);

  // testing

  NumericVector result = leveneBase(x_splits, type);

  double redFlag = result[0];
  result = result[Range(1, 6)];

  if((redFlag == 1.0) & crash) stop("Computation failed");

  if(redFlag == 1.0) return result;

  result.names() =
    CharacterVector::create("n", "k", "f", "df1", "df2", "p_value");

  return result;

}

// version for matrix and a splitting factor

// [[Rcpp::export]]

NumericMatrix leveneMtx(NumericMatrix x,
                        IntegerVector f,
                        String type = "standard",
                        bool crash = true) {

  int nVars = x.ncol();

  NumericMatrix result(nVars, 6);

  colnames(result) =
    CharacterVector::create("n", "k", "f", "df1", "df2", "p_value");

  NumericVector tstOutput(6);

  for(int i = 0; i < nVars; ++i) {

    tstOutput = leveneVec(x(_, i), f, type, crash);

    result(i, _) = tstOutput;

  }

  if(checkNames(x)[1]) rownames(result) = colnames(x);

  return result;

}

// END
