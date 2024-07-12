#ifndef __covariance__
#define __covariance__

#include <Rcpp.h>

using namespace Rcpp;

// Pearson's and Spearman's covariance coefficients for vectors

NumericVector Cov(NumericVector x,
                  NumericVector y,
                  String method);

// concordant and discordant pairs for tau coefficients

NumericVector countPairs(NumericVector x, NumericVector y);

// correlation coefficients for vectors

double rRho(NumericVector x, NumericVector y);
double tauA(NumericVector x, NumericVector y);
double tauB(NumericVector x, NumericVector y);
double xiA(NumericVector x, NumericVector y);
double xiB(NumericVector x, NumericVector y);

NumericVector Cor(NumericVector x,
                  NumericVector y,
                  String method);

#endif // __covariance__
