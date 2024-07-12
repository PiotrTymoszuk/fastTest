/*** Utilities for handling sets */

#include <Rcpp.h>
#include <Rmath.h>
#include <random>


using namespace Rcpp;

// counting elements in intersection and union

NumericVector interSize(CharacterVector x, List dict) {

  int n = dict.size();

  NumericVector counts(n);

  for(int i = 0; i < n; ++i) {

    CharacterVector record = dict[i];

    counts[i] = intersect(x, record).size() * 1.0;

  }

 return counts;

}

NumericVector unionSize(CharacterVector x, List dict) {

  int n = dict.size();

  NumericVector counts(n);

  for(int i = 0; i < n; ++i) {

    CharacterVector record = dict[i];

    counts[i] = union_(x, record).size() * 1.0;

  }

  return counts;

}

// random samples of a character vector

CharacterVector sampleChrVector(CharacterVector x,
                                int size,
                                bool replace = false) {

  int n = x.size();

  if(size > n) stop("Size of the sample cannot exceed the vector size");

  IntegerVector idx = seq(0, n - 1);

  IntegerVector sample_idx = sample(idx, size, replace);

  return x[sample_idx];

}

// a custom 2x2 contingency tables for dictionary entries

// [[Rcpp::export]]

NumericMatrix setCtg(CharacterVector x,
                     CharacterVector all,
                     CharacterVector entry) {

  x = unique(x);
  all = unique(all);
  entry = unique(entry);

  x = na_omit(x);
  all = na_omit(all);
  entry = na_omit(entry);

  double x_size = x.size() * 1.0;
  double all_size = all.size() * 1.0;

  double x_entry_size = intersect(x, entry).size() * 1.0;
  double all_entry_size = intersect(all, entry).size() * 1.0;

  NumericMatrix ctg(2, 2);

  ctg(0, 0) = x_entry_size;
  ctg(1, 0) = all_entry_size;
  ctg(0, 1) = x_size - x_entry_size;
  ctg(1, 1) = all_size - all_entry_size;

  colnames(ctg) = CharacterVector::create("entry", "non_entry");
  rownames(ctg) = CharacterVector::create("x", "all");

  return ctg;

}


// END
