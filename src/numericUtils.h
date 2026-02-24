#ifndef __numericUtils__
#define __numericUtils__

#include <Rcpp.h>

using namespace Rcpp;

// common descriptive stats not implemented by Rcpp

double Median(NumericVector x, bool na_rm);
double Var(NumericVector x, bool na_rm);
double SD(NumericVector x, bool na_rm);
double SEM(NumericVector x, bool na_rm);

// pooled variances used for calculation of Cohen's d

double poolVarPaired(NumericVector x, NumericVector y);
double poolVarWelch(NumericVector x, NumericVector y);
double poolVarStandard(NumericVector x, NumericVector y);

// Hodges - Lehman estimates of location

double locationStd(NumericVector x, NumericVector y);
double locationPair(NumericVector x);

// confidence intervals

NumericVector Quantile(NumericVector x, NumericVector probs);

NumericVector perCI(NumericVector theta, double conf_level);
NumericVector BCA(NumericVector theta, double conf_level);

// Smirnov distribution

NumericVector pSmirnov(double q, double n1, double n2, String alternative);

#endif // __numericUtils__
