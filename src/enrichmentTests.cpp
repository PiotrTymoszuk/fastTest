/*** Set enrichment testing functions */

#include <Rcpp.h>
#include <Rmath.h>
#include <random>
#include "numericUtils.h"
#include "vectorUtils.h"
#include "setUtils.h"

using namespace Rcpp;

// testing for set enrichment with Fisher's exact test:

// Fisher's exact test for a simple 2x2 table

NumericVector setFisherTst(NumericMatrix table,
                           double laplace = 0.0) {

  double x_entry = table(0, 0); // intersect of x and entry
  double all_entry = table(1, 0); // intersect of all and entry
  double x_non_entry = table(0, 1); // x features not in the entry
  double all_non_entry = table(1, 1); // all features not in the entry

  // Calculate the odds ratio

  double odds_ratio =
    ((x_entry + laplace)/x_non_entry)/((all_entry + laplace)/all_non_entry);

  // Calculate the hypergeometric probabilities, twp-tailed case

  double p_value =
    R::phyper(x_entry,
              all_entry,
              all_non_entry,
              x_entry + x_non_entry,
              false, false);

  double p2 =
    R::phyper(x_entry,
              all_entry,
              all_non_entry,
              x_entry + x_non_entry,
              true, false);

  if(p2 < p_value) p_value = p2;

  p_value = 2 * p_value;

  // output

  NumericVector result =
    NumericVector::create(all_entry,
                          x_entry + x_non_entry,
                          x_entry,
                          odds_ratio,
                          p_value);

  result.names() =
    CharacterVector::create("n_entry",
                            "n_x_total",
                            "n_intersect",
                            "or",
                            "p_value");

  return result;

}

// the master function for serial Fisher's test

// [[Rcpp::export]]

NumericMatrix setFisher(CharacterVector x,
                        CharacterVector all,
                        List dict,
                        double laplace = 1.0) {

  // entry control and the result container

  x = na_omit(x);
  all = na_omit(all);

  int dictSize = dict.size();

  NumericMatrix result(dictSize, 5);

  colnames(result) =
    CharacterVector::create("n_entry",
                            "n_x_total",
                            "n_intersect",
                            "or",
                            "p_value");

  // serial testing

  NumericMatrix ctg(2, 2);

  for(int i = 0; i < dictSize; ++i) {

    CharacterVector entry = dict[i];

    ctg = setCtg(x, all, entry);

    result(i, _) = setFisherTst(ctg, laplace);

  }

  return result;

}

// random sample stats

NumericVector getRandStats(const double xCount,
                           NumericVector randCounts,
                           const double xSize,
                           const double entrySize,
                           String ci_type = "bca",
                           double conf_level = 0.95,
                           double laplace = 1.0) {

  // entry control and containers

  randCounts = na_omit(randCounts);
  int randSize = randCounts.size();

  // calculation of OR with confidence intervals

  NumericVector sampleOR = (xCount + laplace)/(randCounts + laplace);

  double meanOR = mean(sampleOR);

  NumericVector CI(2);

  if(ci_type == "bca") {

    CI = BCA(sampleOR, conf_level);

  } else {

    CI = perCI(sampleOR, conf_level);

  }

  // events in favor of the H0 and H1 hypotheses

  NumericVector h0_vector(randSize);

  if(meanOR > 1) {

    h0_vector = ifelse(sampleOR <= 1, 1.0, 0.0);

  } else {

    h0_vector = ifelse(sampleOR >= 1, 1.0, 0.0);

  }

  double h0_number = sum(h0_vector);

  double h1_number = randSize - h0_number;

  double p_value;

  if(h0_number == 0) {

    p_value = (h0_number + 1)/randSize;

  } else {

    p_value = h0_number/randSize;

  }

  // the output vector

  NumericVector outVec =
    NumericVector::create(entrySize, xSize, xCount,
                          meanOR, CI[0], CI[1], randSize,
                          h0_number, h1_number, p_value);

  return outVec;

}


// set enrichment investigated by comparing frequencies in the sample x
// with random samples from the universe vector

// [[Rcpp::export]]

NumericMatrix setEnrichment(CharacterVector x,
                            CharacterVector all,
                            List dict,
                            String ci_type = "bca",
                            double conf_level = 0.95,
                            int n_iter = 100,
                            double laplace = 1.0) {

  // containers

  int x_size = x.size();
  int dict_size = dict.size();

  NumericMatrix results(dict_size, 10);

  colnames(results) =
    CharacterVector::create("n_entry", "n_x_total", "n_intersect",
                            "or", "lower_ci", "upper_ci",
                            "iter_number", "h0_number", "h1_number",
                            "p_value");

  // checking intersects of the input vector with the dictionary list element;
  // sizes of the dictionary entries that overlap with all available features

  NumericVector x_counts = interSize(x, dict);

  NumericVector entry_sizes(dict_size);

  for(int i = 0; i < dict_size; ++i) {

    CharacterVector dict_entry = dict[i];

    entry_sizes[i] = 1.0 * intersect(dict_entry, all).size();

  }

  // random samples of size of the vector x drawn from the all vector
  // checking their intersections with the dictionary elements

  NumericMatrix rand_counts(n_iter, dict_size);
  CharacterVector allSample(x_size);

  for(int i = 0; i < n_iter; ++i) {

    allSample = sampleChrVector(all, x_size, false);

    rand_counts(i, _) = interSize(allSample, dict);

  }

  // counting the cases in support of the hypotheses, OR and CI calculation

  for(int i = 0; i < dict_size; ++i) {

   results(i, _) =
     getRandStats(x_counts[i],
                  rand_counts(_, i),
                  x_size,
                  entry_sizes[i],
                  ci_type,
                  conf_level,
                  laplace);

  }

  return results;

}

// END
