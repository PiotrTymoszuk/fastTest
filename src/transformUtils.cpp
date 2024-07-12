/*** Helper functions for transformation of vectors and matrices */

#include <Rcpp.h>
#include <Rmath.h>
#include "contingencyUtils.h"
#include "transformUtils.h"

using namespace Rcpp;

// checking if a matrix has column and row names

LogicalVector checkNames(const NumericMatrix &x) {

  List s = x.attr("dimnames");  // could be nil or list

  LogicalVector res(2, false);

  if(s.length() == 0) return res;

  if(!Rf_isNull(s[0])) res[0] = true;
  if(!Rf_isNull(s[1])) res[1] = true;

  return res;

}

// checking if a list has names

bool checkListNames(const List &x) {

  if(Rf_isNull(x.attr("names"))) return false;

  return true;

}

// complete cases for a vector pair

List completeCases(NumericVector x, IntegerVector f) {

  int fLen = f.length();

  if(x.length() != fLen) stop("Improper length of f");

  LogicalVector indexes(fLen, true);

  indexes = !(is_na(x) | is_na(f));

  List res(2);

  res[0] = x[indexes];
  res[1] = f[indexes];

  return res;

}

// fast splitting of a vector or a matrix by a vector
// serving as a splitting factor. Positions of f with NAs are omitted from
// from the output

// [[Rcpp::export]]

List Split(NumericVector x, IntegerVector f) {

  // splitting

  IntegerVector f_cat = unique(f).sort();
  int fLen = f.length();

  List res(f_cat.length());

  if(f_cat.length() == 1) {

    res[0] = x;
    res.names() = f_cat[0];

    return res;

  }

  if(x.length() != fLen) stop("Improper length of f");

  for(int k = 0; k < f_cat.length(); ++k) {

    LogicalVector indexes(fLen, false);

    for(int i = 0; i < fLen; ++i) {

      if(f[i] == f_cat[k]) indexes[i] = true;

    }

    // omitting the NA-storing positions from the output

    LogicalVector na_check = !is_na(f);

    indexes = na_check & indexes;

    res[k] = x[indexes];

  }

  res.names() = f_cat;

  return res;

}

// [[Rcpp::export]]

List SplitMtx(NumericMatrix x, IntegerVector f) {

  IntegerVector f_cat = unique(f).sort();

  List res(f_cat.length());

  // a special case handling

  if(f_cat.length() == 1) {

    res[0] = x;
    res.names() = f_cat[0];

    return res;

  }

  // retrieval of row indexes and removal of NA-storing positions

  if(f.length() != x.nrow()) stop("Improper length of f");

  List index_lst(f_cat.length());

  LogicalVector na_check = !is_na(f);

  for(int k = 0; k < f_cat.length(); ++k) {

    LogicalVector indexes_lgl(f.length(), false);
    IntegerVector indexes_int = seq(0, f.length() - 1);

    for(int i = 0; i < f.length(); ++i) {

      if(f[i] == f_cat[k]) indexes_lgl[i] = true;

    }

    indexes_lgl = na_check & indexes_lgl;

    index_lst[k] = indexes_int[indexes_lgl];

  }

  // sub-setting of the matrix

  for(int k = 0; k < f_cat.length(); ++k) {

    IntegerVector row_indexes = index_lst[k];
    NumericMatrix mtx_part(row_indexes.length(), x.ncol());

    for(int i = 0; i < row_indexes.length(); ++i) {

      mtx_part(i, _) = x(i, _);

    }

    res[k] = mtx_part;

    if(checkNames(x)[1]) {

      CharacterVector new_cols = colnames(x);

      colnames(res[k]) = new_cols;

    }

    if(checkNames(x)[0]) {

      CharacterVector new_names = rownames(x);

      rownames(res[k]) = new_names[row_indexes];

    }

  }

  res.names() = f_cat;

  return res;

}

// concatenation of two numeric vectors

NumericVector Concat(NumericVector x, NumericVector y) {

  int xLen = x.length();
  int yLen = y.length();

  NumericVector newVec(xLen + yLen);

  int i = 0;
  int j = 0;

  while(i < xLen + yLen) {

    while(j < xLen) {

      newVec[i] = x[j];

      i += 1;
      j += 1;

    }

    j = 0;

    while(j < yLen) {

      newVec[i] = y[j];

      i += 1;
      j += 1;

    }

  }

  return newVec;

}

// vector pre-processing and checks: length compatibility, non-NA observations
// used by T tests and Wilcoxon tests

List procVectors(NumericVector x,
                 NumericVector y,
                 bool paired = false) {

  // result container

  bool redFlags = false;

  int xLen = x.length();
  int yLen = y.length();

  if(paired) {

    if(xLen == yLen) {

      LogicalVector indexes = !(is_na(x) | is_na(y));

      x = x[indexes];
      y = y[indexes];

    } else {

      x = na_omit(x);
      y = na_omit(y);

      Rcout << "Incompatible length of the input vectors" << "\n";

      redFlags = true;

    }

  } else {

    x = na_omit(x);
    y = na_omit(y);

  }

  xLen = x.length();
  yLen = y.length();

  if((xLen < 2) | (yLen < 2)) {

    Rcout << "Not enough observations to evaluate" << "\n";

    redFlags = true;

  }

  return List::create(x, y, redFlags);


}

// check if categories of a splitting vector are of equal lengths

bool checkEqualSplits(NumericVector x) {

  if(x.length() == 1) return true;

  IntegerVector x_tab = Table(x);

  for(int i = 1; i < x_tab.length(); ++i) {

    if(x_tab[i] != x_tab[i - 1]) return false;

  }

  return true;

}

// END
