/*** General calculation helpers used by other Cpp functions */

#include <Rcpp.h>
#include <Rmath.h>
#include "rankUtils.h"
#include "numericUtils.h"
#include "vectorUtils.h"
#include "contingencyUtils.h"
#include "transformUtils.h"

using namespace Rcpp;


// Define a custom comparator class to compare values in the vector
// and handling of NAs

class Comparator {

private:

  // definition of the is_na() method for the Comparator class

  const NumericVector& ref;

  bool is_na(double x) const {

    return traits::is_na<REALSXP>(x);

  }

public:

  Comparator(const Rcpp::NumericVector& ref_) : ref(ref_) {}

  // defines comparison operation for the vector elements;
  // NA should land at the end of the ranks!!!

  bool operator()(const int ilhs, const int irhs) const {

    double lhs = ref[ilhs], rhs = ref[irhs];

    if (is_na(lhs)) return false;

    if (is_na(rhs)) return true;

    return lhs < rhs;

  }

};

// ranks of a numeric vector

NumericVector Rank(NumericVector x,
                   String ties_method = "average",
                   bool na_rm = true) {

  if(na_rm) x = na_omit(x);

  // Get the size of the input vector

  R_xlen_t sz = x.length();

  // Create an integer vector to store the sorted indices

  IntegerVector w = seq(0, sz - 1);

  std::sort(w.begin(), w.end(), Comparator(x));

  // Initialize a numeric vector to store the ranks

  NumericVector r = no_init_vector(sz);

  // Iterate through the sorted indices to compute ranks

  for (R_xlen_t n, i = 0; i < sz; i += n) {

    n = 1;

    while (i + n < sz && x[w[i]] == x[w[i + n]]) ++n;

    if(ties_method == "random") {

      IntegerVector freeRanks = seq(i + 1, i + n);

      IntegerVector randRanks = sample(freeRanks, freeRanks.length(), false);

      for (R_xlen_t k = 0; k < n; k++) r[w[i + k]] = randRanks[k];

    } else {

      for (R_xlen_t k = 0; k < n; k++) {

        if (ties_method == "average") {

          r[w[i + k]] = i + (n + 1) / 2.0;  // Average ranks

        } else if(ties_method == "min") {

          r[w[i + k]] = i + 1;  // Minimum ranks

        } else if(ties_method == "max") {

          r[w[i + k]] = i + n; // Maximum ranks

        }

      }

    }

  }

  return r;

}

// signed ranks: zeroes are silently removed

NumericVector signedRank(NumericVector x,
                         String ties_method = "average",
                         bool na_rm = true) {

  // NA and zero removal

  if(na_rm) x = na_omit(x);

  int xLen = x.length();

  LogicalVector non_zero(xLen, true);

  for(int i = 0; i < xLen; ++i) {

    if(x[i] == 0) non_zero[i] = false;

  }

  x = x[non_zero];

  // rank calculation

  NumericVector x_abs = Abs(x, false);
  NumericVector x_sign = Sign(x, false);

  NumericVector r = Rank(x_abs, ties_method, false);

  return product(r, x_sign);

}

// processing of vectors for Kruskal-Wallis and similar tests
// should return a vector of ranks, the splitting factor, and
// a logical value 'redFlag' is returned if there are critical vector
// issues

List rankVectors(NumericVector x, IntegerVector f) {

  // result container

  List result(3);

  result.names() = CharacterVector::create("ranks", "f", "redFlag");

  // complete cases

  LogicalVector idx = !((is_na(x)) | (is_na(f)));

  x = x[idx];
  f = f[idx];

  // ranks

  NumericVector r = Rank(x, "average", false);

  result[0] = r;
  result[1] = f;
  result[2] = false; // by default, there are no issues

  // control of the input: there are basically three crash possibilities:

  if(unique(f).length() < 2) {

    warning("At least two analysis groups are required");

    result[2] = true;

    return result;

  }

  IntegerVector f_tab = Table(as<NumericVector>(f));

  for(int i = 0; i < f_tab.length(); ++i) {

    if(f_tab[i] < 2) {

      warning("Not enough observations in one of the groups");

      result[2] = true;

      return result;

    }

  }

  if(unique(r).length() < 2) {

    warning("All observations are tied");

    result[2] = true;

    return result;

  }

  return result;

}

// processing the input vectors for Friedman test.
// The function returns a list of:
// 1) a numeric matrix of block-wise ranks;
// blocks define the rows, treatments define the columns
// 2) a logical value indicating data problems

// [[Rcpp::export]]

List rankBlockVectors(NumericVector x,
                      IntegerVector f,
                      IntegerVector b) {

  // result container

  List result(2);

  result.names() = CharacterVector::create("rankMtx", "redFlag");

  result[1] = false;

  // complete cases and universal variables

  LogicalVector idx = !((is_na(x)) | (is_na(f)) | (is_na(b)));

  x = x[idx];
  f = f[idx];
  b = b[idx];

  int k = unique(f).length(); // number of treatment groups
  int B = unique(b).length(); // number of blocks

  // checking for possible problems with the data

  bool fEqual = checkEqualSplits(as<NumericVector>(f));
  bool bEqual = checkEqualSplits(as<NumericVector>(b));

  if(!fEqual | !bEqual) {

    warning("Level-deficient design");

    result[1] = true;

    return result;

  }

  if(k < 2) {

    warning("At least two analysis groups are required");

    result[1] = true;

    return result;

  }

  IntegerVector f_tab = Table(as<NumericVector>(f));

  for(int i = 0; i < f_tab.length(); ++i) {

    if(f_tab[i] < 2) {

      warning("Not enough observations in one of the groups");

      result[1] = true;

      return result;

    }

  }

  if(unique(b).length() < 2) {

    warning("At least two analysis groups are required");

    result[1] = true;

    return result;

  }


  if(B < 2) {

    warning("At least two blocks are required");

    result[1] = true;

    return result;

  }

  IntegerVector b_tab = Table(as<NumericVector>(b));

  for(int i = 0; i < b_tab.length(); ++i) {

    if(f_tab[i] < 2) {

      warning("Not enough observations in one of the blocks");

      result[3] = true;

      return result;

    }

  }

  // ranks within the blocks

  List x_splits = Split(x, b);

  NumericMatrix rankMtx(B, k);

  for(int i = 0; i < B; ++i) {

    NumericVector groupRanks = Rank(x_splits[i], "average", false);

    rankMtx(i, _) = groupRanks;

  }

  result[0] = rankMtx;

  return result;

}

// END
