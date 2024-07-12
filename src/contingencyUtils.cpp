/*** Contingency tables */

#include <Rcpp.h>
#include <Rmath.h>
#include "contingencyUtils.h"
#include "vectorUtils.h"
#include "transformUtils.h"

using namespace Rcpp;

// calculation of element frequency in a numeric vector
// NAs are skipped!
// Table creates a an equivalent of R's table (i.e. a named vector of
// frequencies), while xTable produces a contingency matrix

// [[Rcpp::export]]

IntegerVector Table(NumericVector x) {

  IntegerVector counts;

  x = na_omit(x);

  counts = table(x);

  return counts;

}

// [[Rcpp::export]]

NumericMatrix xTable(NumericVector x, IntegerVector f) {

  // checking for complete cases

   List compCases = completeCases(x, f);

    x = compCases[0];
    f = compCases[1];

  // vectors for categories of x and f and their counts

  NumericVector x_cat = unique(x).sort();
  IntegerVector f_cat = unique(f).sort();

  int x_dim = x_cat.length();
  int f_dim = f_cat.length();

  // vectors storing the subsetting indexes and tallying results

  List x_splits(f_dim);

  // a result matrix with all available categories of x in columns
  // and categories of f in rows

  NumericMatrix res(x_dim, f_dim);

  // a special case handling

  if(f_dim == 1) {

    IntegerVector res_col = Table(x);
    res(_, 0) = res_col;

    return res;

  }

  // typical case: multiple categories of f

  if(x.length() != f.length()) stop("Improper length of f");

  x_splits = Split(x, f);

  for(int j = 0; j < f_dim; ++j) {

    int counter = 0;
    NumericVector x_part;

    for(int i = 0; i < x_dim; ++i) {

      x_part = x_splits[j];

      counter = sum(ifelse(x_part == x_cat[i], 1, 0));

      res(i, j) = counter;

    }

  }

  return res;

}

// END
