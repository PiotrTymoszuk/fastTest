// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// chiSqTestTbl
NumericVector chiSqTestTbl(NumericMatrix ctg, bool correct, bool crash);
RcppExport SEXP _fastTest_chiSqTestTbl(SEXP ctgSEXP, SEXP correctSEXP, SEXP crashSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ctg(ctgSEXP);
    Rcpp::traits::input_parameter< bool >::type correct(correctSEXP);
    Rcpp::traits::input_parameter< bool >::type crash(crashSEXP);
    rcpp_result_gen = Rcpp::wrap(chiSqTestTbl(ctg, correct, crash));
    return rcpp_result_gen;
END_RCPP
}
// chiSqTestVec
NumericVector chiSqTestVec(NumericVector x, IntegerVector f, bool correct, bool crash);
RcppExport SEXP _fastTest_chiSqTestVec(SEXP xSEXP, SEXP fSEXP, SEXP correctSEXP, SEXP crashSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< bool >::type correct(correctSEXP);
    Rcpp::traits::input_parameter< bool >::type crash(crashSEXP);
    rcpp_result_gen = Rcpp::wrap(chiSqTestVec(x, f, correct, crash));
    return rcpp_result_gen;
END_RCPP
}
// chiSqTestMtx
NumericMatrix chiSqTestMtx(NumericMatrix x, IntegerVector f, bool correct, bool crash);
RcppExport SEXP _fastTest_chiSqTestMtx(SEXP xSEXP, SEXP fSEXP, SEXP correctSEXP, SEXP crashSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< bool >::type correct(correctSEXP);
    Rcpp::traits::input_parameter< bool >::type crash(crashSEXP);
    rcpp_result_gen = Rcpp::wrap(chiSqTestMtx(x, f, correct, crash));
    return rcpp_result_gen;
END_RCPP
}
// Table
IntegerVector Table(NumericVector x);
RcppExport SEXP _fastTest_Table(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(Table(x));
    return rcpp_result_gen;
END_RCPP
}
// xTable
NumericMatrix xTable(NumericVector x, IntegerVector f);
RcppExport SEXP _fastTest_xTable(SEXP xSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(xTable(x, f));
    return rcpp_result_gen;
END_RCPP
}
// permCorVec
NumericVector permCorVec(NumericVector x, NumericVector y, String method, String alternative, int n_iter);
RcppExport SEXP _fastTest_permCorVec(SEXP xSEXP, SEXP ySEXP, SEXP methodSEXP, SEXP alternativeSEXP, SEXP n_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< String >::type method(methodSEXP);
    Rcpp::traits::input_parameter< String >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< int >::type n_iter(n_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(permCorVec(x, y, method, alternative, n_iter));
    return rcpp_result_gen;
END_RCPP
}
// permCorMtx
NumericMatrix permCorMtx(NumericMatrix x, String method, String alternative, int n_iter);
RcppExport SEXP _fastTest_permCorMtx(SEXP xSEXP, SEXP methodSEXP, SEXP alternativeSEXP, SEXP n_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< String >::type method(methodSEXP);
    Rcpp::traits::input_parameter< String >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< int >::type n_iter(n_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(permCorMtx(x, method, alternative, n_iter));
    return rcpp_result_gen;
END_RCPP
}
// bootCorVec
NumericVector bootCorVec(NumericVector x, NumericVector y, String method, String ci_type, double conf_level, int n_iter);
RcppExport SEXP _fastTest_bootCorVec(SEXP xSEXP, SEXP ySEXP, SEXP methodSEXP, SEXP ci_typeSEXP, SEXP conf_levelSEXP, SEXP n_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< String >::type method(methodSEXP);
    Rcpp::traits::input_parameter< String >::type ci_type(ci_typeSEXP);
    Rcpp::traits::input_parameter< double >::type conf_level(conf_levelSEXP);
    Rcpp::traits::input_parameter< int >::type n_iter(n_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(bootCorVec(x, y, method, ci_type, conf_level, n_iter));
    return rcpp_result_gen;
END_RCPP
}
// bootCorMtx
NumericVector bootCorMtx(NumericMatrix x, String method, String ci_type, double conf_level, int n_iter);
RcppExport SEXP _fastTest_bootCorMtx(SEXP xSEXP, SEXP methodSEXP, SEXP ci_typeSEXP, SEXP conf_levelSEXP, SEXP n_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< String >::type method(methodSEXP);
    Rcpp::traits::input_parameter< String >::type ci_type(ci_typeSEXP);
    Rcpp::traits::input_parameter< double >::type conf_level(conf_levelSEXP);
    Rcpp::traits::input_parameter< int >::type n_iter(n_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(bootCorMtx(x, method, ci_type, conf_level, n_iter));
    return rcpp_result_gen;
END_RCPP
}
// Cov
NumericVector Cov(NumericVector x, NumericVector y, String method);
RcppExport SEXP _fastTest_Cov(SEXP xSEXP, SEXP ySEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< String >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(Cov(x, y, method));
    return rcpp_result_gen;
END_RCPP
}
// CovMtx
NumericMatrix CovMtx(NumericMatrix x, String method);
RcppExport SEXP _fastTest_CovMtx(SEXP xSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< String >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(CovMtx(x, method));
    return rcpp_result_gen;
END_RCPP
}
// Cor
NumericVector Cor(NumericVector x, NumericVector y, String method);
RcppExport SEXP _fastTest_Cor(SEXP xSEXP, SEXP ySEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< String >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(Cor(x, y, method));
    return rcpp_result_gen;
END_RCPP
}
// CorMtx
NumericMatrix CorMtx(NumericMatrix x, String method);
RcppExport SEXP _fastTest_CorMtx(SEXP xSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< String >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(CorMtx(x, method));
    return rcpp_result_gen;
END_RCPP
}
// setFisher
NumericMatrix setFisher(CharacterVector x, CharacterVector all, List dict, double laplace);
RcppExport SEXP _fastTest_setFisher(SEXP xSEXP, SEXP allSEXP, SEXP dictSEXP, SEXP laplaceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type all(allSEXP);
    Rcpp::traits::input_parameter< List >::type dict(dictSEXP);
    Rcpp::traits::input_parameter< double >::type laplace(laplaceSEXP);
    rcpp_result_gen = Rcpp::wrap(setFisher(x, all, dict, laplace));
    return rcpp_result_gen;
END_RCPP
}
// setEnrichment
NumericMatrix setEnrichment(CharacterVector x, CharacterVector all, List dict, String ci_type, double conf_level, int n_iter, double laplace);
RcppExport SEXP _fastTest_setEnrichment(SEXP xSEXP, SEXP allSEXP, SEXP dictSEXP, SEXP ci_typeSEXP, SEXP conf_levelSEXP, SEXP n_iterSEXP, SEXP laplaceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type all(allSEXP);
    Rcpp::traits::input_parameter< List >::type dict(dictSEXP);
    Rcpp::traits::input_parameter< String >::type ci_type(ci_typeSEXP);
    Rcpp::traits::input_parameter< double >::type conf_level(conf_levelSEXP);
    Rcpp::traits::input_parameter< int >::type n_iter(n_iterSEXP);
    Rcpp::traits::input_parameter< double >::type laplace(laplaceSEXP);
    rcpp_result_gen = Rcpp::wrap(setEnrichment(x, all, dict, ci_type, conf_level, n_iter, laplace));
    return rcpp_result_gen;
END_RCPP
}
// friedmanVec
NumericVector friedmanVec(NumericVector x, IntegerVector f, IntegerVector b, bool crash);
RcppExport SEXP _fastTest_friedmanVec(SEXP xSEXP, SEXP fSEXP, SEXP bSEXP, SEXP crashSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< bool >::type crash(crashSEXP);
    rcpp_result_gen = Rcpp::wrap(friedmanVec(x, f, b, crash));
    return rcpp_result_gen;
END_RCPP
}
// friedmanMtx
NumericMatrix friedmanMtx(NumericMatrix x, IntegerVector f, IntegerVector b, bool crash);
RcppExport SEXP _fastTest_friedmanMtx(SEXP xSEXP, SEXP fSEXP, SEXP bSEXP, SEXP crashSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< bool >::type crash(crashSEXP);
    rcpp_result_gen = Rcpp::wrap(friedmanMtx(x, f, b, crash));
    return rcpp_result_gen;
END_RCPP
}
// kruskalVec
NumericVector kruskalVec(NumericVector x, IntegerVector f, bool crash);
RcppExport SEXP _fastTest_kruskalVec(SEXP xSEXP, SEXP fSEXP, SEXP crashSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< bool >::type crash(crashSEXP);
    rcpp_result_gen = Rcpp::wrap(kruskalVec(x, f, crash));
    return rcpp_result_gen;
END_RCPP
}
// kruskalMtx
NumericVector kruskalMtx(NumericMatrix x, IntegerVector f, bool crash);
RcppExport SEXP _fastTest_kruskalMtx(SEXP xSEXP, SEXP fSEXP, SEXP crashSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< bool >::type crash(crashSEXP);
    rcpp_result_gen = Rcpp::wrap(kruskalMtx(x, f, crash));
    return rcpp_result_gen;
END_RCPP
}
// leveneBase
NumericVector leveneBase(List x, String type);
RcppExport SEXP _fastTest_leveneBase(SEXP xSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< String >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(leveneBase(x, type));
    return rcpp_result_gen;
END_RCPP
}
// leveneVec
NumericVector leveneVec(NumericVector x, IntegerVector f, String type, bool crash);
RcppExport SEXP _fastTest_leveneVec(SEXP xSEXP, SEXP fSEXP, SEXP typeSEXP, SEXP crashSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< String >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type crash(crashSEXP);
    rcpp_result_gen = Rcpp::wrap(leveneVec(x, f, type, crash));
    return rcpp_result_gen;
END_RCPP
}
// leveneMtx
NumericMatrix leveneMtx(NumericMatrix x, IntegerVector f, String type, bool crash);
RcppExport SEXP _fastTest_leveneMtx(SEXP xSEXP, SEXP fSEXP, SEXP typeSEXP, SEXP crashSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< String >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type crash(crashSEXP);
    rcpp_result_gen = Rcpp::wrap(leveneMtx(x, f, type, crash));
    return rcpp_result_gen;
END_RCPP
}
// metaVec
NumericVector metaVec(NumericVector y, NumericVector e, String type, String alternative, double conf_level, bool crash);
RcppExport SEXP _fastTest_metaVec(SEXP ySEXP, SEXP eSEXP, SEXP typeSEXP, SEXP alternativeSEXP, SEXP conf_levelSEXP, SEXP crashSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type e(eSEXP);
    Rcpp::traits::input_parameter< String >::type type(typeSEXP);
    Rcpp::traits::input_parameter< String >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< double >::type conf_level(conf_levelSEXP);
    Rcpp::traits::input_parameter< bool >::type crash(crashSEXP);
    rcpp_result_gen = Rcpp::wrap(metaVec(y, e, type, alternative, conf_level, crash));
    return rcpp_result_gen;
END_RCPP
}
// metaMtx
NumericMatrix metaMtx(NumericMatrix y, NumericMatrix e, String type, String alternative, double conf_level, bool crash);
RcppExport SEXP _fastTest_metaMtx(SEXP ySEXP, SEXP eSEXP, SEXP typeSEXP, SEXP alternativeSEXP, SEXP conf_levelSEXP, SEXP crashSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type e(eSEXP);
    Rcpp::traits::input_parameter< String >::type type(typeSEXP);
    Rcpp::traits::input_parameter< String >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< double >::type conf_level(conf_levelSEXP);
    Rcpp::traits::input_parameter< bool >::type crash(crashSEXP);
    rcpp_result_gen = Rcpp::wrap(metaMtx(y, e, type, alternative, conf_level, crash));
    return rcpp_result_gen;
END_RCPP
}
// metaList
NumericMatrix metaList(List y, List e, String type, String alternative, double conf_level, bool crash);
RcppExport SEXP _fastTest_metaList(SEXP ySEXP, SEXP eSEXP, SEXP typeSEXP, SEXP alternativeSEXP, SEXP conf_levelSEXP, SEXP crashSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type y(ySEXP);
    Rcpp::traits::input_parameter< List >::type e(eSEXP);
    Rcpp::traits::input_parameter< String >::type type(typeSEXP);
    Rcpp::traits::input_parameter< String >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< double >::type conf_level(conf_levelSEXP);
    Rcpp::traits::input_parameter< bool >::type crash(crashSEXP);
    rcpp_result_gen = Rcpp::wrap(metaList(y, e, type, alternative, conf_level, crash));
    return rcpp_result_gen;
END_RCPP
}
// Median
double Median(NumericVector x, bool na_rm);
RcppExport SEXP _fastTest_Median(SEXP xSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(Median(x, na_rm));
    return rcpp_result_gen;
END_RCPP
}
// Var
double Var(NumericVector x, bool na_rm);
RcppExport SEXP _fastTest_Var(SEXP xSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(Var(x, na_rm));
    return rcpp_result_gen;
END_RCPP
}
// SD
double SD(NumericVector x, bool na_rm);
RcppExport SEXP _fastTest_SD(SEXP xSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(SD(x, na_rm));
    return rcpp_result_gen;
END_RCPP
}
// SEM
double SEM(NumericVector x, bool na_rm);
RcppExport SEXP _fastTest_SEM(SEXP xSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(SEM(x, na_rm));
    return rcpp_result_gen;
END_RCPP
}
// locationStd
double locationStd(NumericVector x, NumericVector y);
RcppExport SEXP _fastTest_locationStd(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(locationStd(x, y));
    return rcpp_result_gen;
END_RCPP
}
// locationPair
double locationPair(NumericVector x);
RcppExport SEXP _fastTest_locationPair(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(locationPair(x));
    return rcpp_result_gen;
END_RCPP
}
// Quantile
NumericVector Quantile(NumericVector x, NumericVector probs);
RcppExport SEXP _fastTest_Quantile(SEXP xSEXP, SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type probs(probsSEXP);
    rcpp_result_gen = Rcpp::wrap(Quantile(x, probs));
    return rcpp_result_gen;
END_RCPP
}
// perCI
NumericVector perCI(NumericVector theta, double conf_level);
RcppExport SEXP _fastTest_perCI(SEXP thetaSEXP, SEXP conf_levelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type conf_level(conf_levelSEXP);
    rcpp_result_gen = Rcpp::wrap(perCI(theta, conf_level));
    return rcpp_result_gen;
END_RCPP
}
// BCA
NumericVector BCA(NumericVector theta, double conf_level);
RcppExport SEXP _fastTest_BCA(SEXP thetaSEXP, SEXP conf_levelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type conf_level(conf_levelSEXP);
    rcpp_result_gen = Rcpp::wrap(BCA(theta, conf_level));
    return rcpp_result_gen;
END_RCPP
}
// oneAnovaBase
NumericVector oneAnovaBase(List x);
RcppExport SEXP _fastTest_oneAnovaBase(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(oneAnovaBase(x));
    return rcpp_result_gen;
END_RCPP
}
// oneAnovaVec
NumericVector oneAnovaVec(NumericVector x, IntegerVector f, bool crash);
RcppExport SEXP _fastTest_oneAnovaVec(SEXP xSEXP, SEXP fSEXP, SEXP crashSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< bool >::type crash(crashSEXP);
    rcpp_result_gen = Rcpp::wrap(oneAnovaVec(x, f, crash));
    return rcpp_result_gen;
END_RCPP
}
// oneAnovaMtx
NumericMatrix oneAnovaMtx(NumericMatrix x, IntegerVector f, bool crash);
RcppExport SEXP _fastTest_oneAnovaMtx(SEXP xSEXP, SEXP fSEXP, SEXP crashSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< bool >::type crash(crashSEXP);
    rcpp_result_gen = Rcpp::wrap(oneAnovaMtx(x, f, crash));
    return rcpp_result_gen;
END_RCPP
}
// oneBlockAnovaVec
NumericVector oneBlockAnovaVec(NumericVector x, IntegerVector f, IntegerVector b);
RcppExport SEXP _fastTest_oneBlockAnovaVec(SEXP xSEXP, SEXP fSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(oneBlockAnovaVec(x, f, b));
    return rcpp_result_gen;
END_RCPP
}
// oneBlockAnovaMtx
NumericMatrix oneBlockAnovaMtx(NumericMatrix x, IntegerVector f, IntegerVector b, bool crash);
RcppExport SEXP _fastTest_oneBlockAnovaMtx(SEXP xSEXP, SEXP fSEXP, SEXP bSEXP, SEXP crashSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< bool >::type crash(crashSEXP);
    rcpp_result_gen = Rcpp::wrap(oneBlockAnovaMtx(x, f, b, crash));
    return rcpp_result_gen;
END_RCPP
}
// rankBlockVectors
List rankBlockVectors(NumericVector x, IntegerVector f, IntegerVector b);
RcppExport SEXP _fastTest_rankBlockVectors(SEXP xSEXP, SEXP fSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(rankBlockVectors(x, f, b));
    return rcpp_result_gen;
END_RCPP
}
// setCtg
NumericMatrix setCtg(CharacterVector x, CharacterVector all, CharacterVector entry);
RcppExport SEXP _fastTest_setCtg(SEXP xSEXP, SEXP allSEXP, SEXP entrySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type all(allSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type entry(entrySEXP);
    rcpp_result_gen = Rcpp::wrap(setCtg(x, all, entry));
    return rcpp_result_gen;
END_RCPP
}
// tTestStd
NumericVector tTestStd(NumericVector x, NumericVector y, String type, String alternative, double conf_level, bool crash);
RcppExport SEXP _fastTest_tTestStd(SEXP xSEXP, SEXP ySEXP, SEXP typeSEXP, SEXP alternativeSEXP, SEXP conf_levelSEXP, SEXP crashSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< String >::type type(typeSEXP);
    Rcpp::traits::input_parameter< String >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< double >::type conf_level(conf_levelSEXP);
    Rcpp::traits::input_parameter< bool >::type crash(crashSEXP);
    rcpp_result_gen = Rcpp::wrap(tTestStd(x, y, type, alternative, conf_level, crash));
    return rcpp_result_gen;
END_RCPP
}
// tTestVec
NumericVector tTestVec(NumericVector x, IntegerVector f, String type, String alternative, double conf_level, bool crash);
RcppExport SEXP _fastTest_tTestVec(SEXP xSEXP, SEXP fSEXP, SEXP typeSEXP, SEXP alternativeSEXP, SEXP conf_levelSEXP, SEXP crashSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< String >::type type(typeSEXP);
    Rcpp::traits::input_parameter< String >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< double >::type conf_level(conf_levelSEXP);
    Rcpp::traits::input_parameter< bool >::type crash(crashSEXP);
    rcpp_result_gen = Rcpp::wrap(tTestVec(x, f, type, alternative, conf_level, crash));
    return rcpp_result_gen;
END_RCPP
}
// tTestMtx
NumericMatrix tTestMtx(NumericMatrix x, IntegerVector f, String type, String alternative, double conf_level, bool crash);
RcppExport SEXP _fastTest_tTestMtx(SEXP xSEXP, SEXP fSEXP, SEXP typeSEXP, SEXP alternativeSEXP, SEXP conf_levelSEXP, SEXP crashSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< String >::type type(typeSEXP);
    Rcpp::traits::input_parameter< String >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< double >::type conf_level(conf_levelSEXP);
    Rcpp::traits::input_parameter< bool >::type crash(crashSEXP);
    rcpp_result_gen = Rcpp::wrap(tTestMtx(x, f, type, alternative, conf_level, crash));
    return rcpp_result_gen;
END_RCPP
}
// Split
List Split(NumericVector x, IntegerVector f);
RcppExport SEXP _fastTest_Split(SEXP xSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(Split(x, f));
    return rcpp_result_gen;
END_RCPP
}
// SplitMtx
List SplitMtx(NumericMatrix x, IntegerVector f);
RcppExport SEXP _fastTest_SplitMtx(SEXP xSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(SplitMtx(x, f));
    return rcpp_result_gen;
END_RCPP
}
// outerSum
NumericMatrix outerSum(NumericVector x, NumericVector y, bool na_rm);
RcppExport SEXP _fastTest_outerSum(SEXP xSEXP, SEXP ySEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(outerSum(x, y, na_rm));
    return rcpp_result_gen;
END_RCPP
}
// outerDelta
NumericMatrix outerDelta(NumericVector x, NumericVector y, bool na_rm);
RcppExport SEXP _fastTest_outerDelta(SEXP xSEXP, SEXP ySEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(outerDelta(x, y, na_rm));
    return rcpp_result_gen;
END_RCPP
}
// outerProduct
NumericMatrix outerProduct(NumericVector x, NumericVector y, bool na_rm);
RcppExport SEXP _fastTest_outerProduct(SEXP xSEXP, SEXP ySEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(outerProduct(x, y, na_rm));
    return rcpp_result_gen;
END_RCPP
}
// testWilcoxonVec
NumericVector testWilcoxonVec(NumericVector x, IntegerVector f, String type, String alternative, bool exact, bool correct, bool hl_estimate, bool crash);
RcppExport SEXP _fastTest_testWilcoxonVec(SEXP xSEXP, SEXP fSEXP, SEXP typeSEXP, SEXP alternativeSEXP, SEXP exactSEXP, SEXP correctSEXP, SEXP hl_estimateSEXP, SEXP crashSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< String >::type type(typeSEXP);
    Rcpp::traits::input_parameter< String >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< bool >::type exact(exactSEXP);
    Rcpp::traits::input_parameter< bool >::type correct(correctSEXP);
    Rcpp::traits::input_parameter< bool >::type hl_estimate(hl_estimateSEXP);
    Rcpp::traits::input_parameter< bool >::type crash(crashSEXP);
    rcpp_result_gen = Rcpp::wrap(testWilcoxonVec(x, f, type, alternative, exact, correct, hl_estimate, crash));
    return rcpp_result_gen;
END_RCPP
}
// testWilcoxonStd
NumericVector testWilcoxonStd(NumericVector x, NumericVector y, String type, String alternative, bool exact, bool correct, bool hl_estimate, bool crash);
RcppExport SEXP _fastTest_testWilcoxonStd(SEXP xSEXP, SEXP ySEXP, SEXP typeSEXP, SEXP alternativeSEXP, SEXP exactSEXP, SEXP correctSEXP, SEXP hl_estimateSEXP, SEXP crashSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< String >::type type(typeSEXP);
    Rcpp::traits::input_parameter< String >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< bool >::type exact(exactSEXP);
    Rcpp::traits::input_parameter< bool >::type correct(correctSEXP);
    Rcpp::traits::input_parameter< bool >::type hl_estimate(hl_estimateSEXP);
    Rcpp::traits::input_parameter< bool >::type crash(crashSEXP);
    rcpp_result_gen = Rcpp::wrap(testWilcoxonStd(x, y, type, alternative, exact, correct, hl_estimate, crash));
    return rcpp_result_gen;
END_RCPP
}
// testWilcoxonMtx
NumericMatrix testWilcoxonMtx(NumericMatrix x, IntegerVector f, String type, String alternative, bool exact, bool correct, bool hl_estimate, bool crash);
RcppExport SEXP _fastTest_testWilcoxonMtx(SEXP xSEXP, SEXP fSEXP, SEXP typeSEXP, SEXP alternativeSEXP, SEXP exactSEXP, SEXP correctSEXP, SEXP hl_estimateSEXP, SEXP crashSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< String >::type type(typeSEXP);
    Rcpp::traits::input_parameter< String >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< bool >::type exact(exactSEXP);
    Rcpp::traits::input_parameter< bool >::type correct(correctSEXP);
    Rcpp::traits::input_parameter< bool >::type hl_estimate(hl_estimateSEXP);
    Rcpp::traits::input_parameter< bool >::type crash(crashSEXP);
    rcpp_result_gen = Rcpp::wrap(testWilcoxonMtx(x, f, type, alternative, exact, correct, hl_estimate, crash));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fastTest_chiSqTestTbl", (DL_FUNC) &_fastTest_chiSqTestTbl, 3},
    {"_fastTest_chiSqTestVec", (DL_FUNC) &_fastTest_chiSqTestVec, 4},
    {"_fastTest_chiSqTestMtx", (DL_FUNC) &_fastTest_chiSqTestMtx, 4},
    {"_fastTest_Table", (DL_FUNC) &_fastTest_Table, 1},
    {"_fastTest_xTable", (DL_FUNC) &_fastTest_xTable, 2},
    {"_fastTest_permCorVec", (DL_FUNC) &_fastTest_permCorVec, 5},
    {"_fastTest_permCorMtx", (DL_FUNC) &_fastTest_permCorMtx, 4},
    {"_fastTest_bootCorVec", (DL_FUNC) &_fastTest_bootCorVec, 6},
    {"_fastTest_bootCorMtx", (DL_FUNC) &_fastTest_bootCorMtx, 5},
    {"_fastTest_Cov", (DL_FUNC) &_fastTest_Cov, 3},
    {"_fastTest_CovMtx", (DL_FUNC) &_fastTest_CovMtx, 2},
    {"_fastTest_Cor", (DL_FUNC) &_fastTest_Cor, 3},
    {"_fastTest_CorMtx", (DL_FUNC) &_fastTest_CorMtx, 2},
    {"_fastTest_setFisher", (DL_FUNC) &_fastTest_setFisher, 4},
    {"_fastTest_setEnrichment", (DL_FUNC) &_fastTest_setEnrichment, 7},
    {"_fastTest_friedmanVec", (DL_FUNC) &_fastTest_friedmanVec, 4},
    {"_fastTest_friedmanMtx", (DL_FUNC) &_fastTest_friedmanMtx, 4},
    {"_fastTest_kruskalVec", (DL_FUNC) &_fastTest_kruskalVec, 3},
    {"_fastTest_kruskalMtx", (DL_FUNC) &_fastTest_kruskalMtx, 3},
    {"_fastTest_leveneBase", (DL_FUNC) &_fastTest_leveneBase, 2},
    {"_fastTest_leveneVec", (DL_FUNC) &_fastTest_leveneVec, 4},
    {"_fastTest_leveneMtx", (DL_FUNC) &_fastTest_leveneMtx, 4},
    {"_fastTest_metaVec", (DL_FUNC) &_fastTest_metaVec, 6},
    {"_fastTest_metaMtx", (DL_FUNC) &_fastTest_metaMtx, 6},
    {"_fastTest_metaList", (DL_FUNC) &_fastTest_metaList, 6},
    {"_fastTest_Median", (DL_FUNC) &_fastTest_Median, 2},
    {"_fastTest_Var", (DL_FUNC) &_fastTest_Var, 2},
    {"_fastTest_SD", (DL_FUNC) &_fastTest_SD, 2},
    {"_fastTest_SEM", (DL_FUNC) &_fastTest_SEM, 2},
    {"_fastTest_locationStd", (DL_FUNC) &_fastTest_locationStd, 2},
    {"_fastTest_locationPair", (DL_FUNC) &_fastTest_locationPair, 1},
    {"_fastTest_Quantile", (DL_FUNC) &_fastTest_Quantile, 2},
    {"_fastTest_perCI", (DL_FUNC) &_fastTest_perCI, 2},
    {"_fastTest_BCA", (DL_FUNC) &_fastTest_BCA, 2},
    {"_fastTest_oneAnovaBase", (DL_FUNC) &_fastTest_oneAnovaBase, 1},
    {"_fastTest_oneAnovaVec", (DL_FUNC) &_fastTest_oneAnovaVec, 3},
    {"_fastTest_oneAnovaMtx", (DL_FUNC) &_fastTest_oneAnovaMtx, 3},
    {"_fastTest_oneBlockAnovaVec", (DL_FUNC) &_fastTest_oneBlockAnovaVec, 3},
    {"_fastTest_oneBlockAnovaMtx", (DL_FUNC) &_fastTest_oneBlockAnovaMtx, 4},
    {"_fastTest_rankBlockVectors", (DL_FUNC) &_fastTest_rankBlockVectors, 3},
    {"_fastTest_setCtg", (DL_FUNC) &_fastTest_setCtg, 3},
    {"_fastTest_tTestStd", (DL_FUNC) &_fastTest_tTestStd, 6},
    {"_fastTest_tTestVec", (DL_FUNC) &_fastTest_tTestVec, 6},
    {"_fastTest_tTestMtx", (DL_FUNC) &_fastTest_tTestMtx, 6},
    {"_fastTest_Split", (DL_FUNC) &_fastTest_Split, 2},
    {"_fastTest_SplitMtx", (DL_FUNC) &_fastTest_SplitMtx, 2},
    {"_fastTest_outerSum", (DL_FUNC) &_fastTest_outerSum, 3},
    {"_fastTest_outerDelta", (DL_FUNC) &_fastTest_outerDelta, 3},
    {"_fastTest_outerProduct", (DL_FUNC) &_fastTest_outerProduct, 3},
    {"_fastTest_testWilcoxonVec", (DL_FUNC) &_fastTest_testWilcoxonVec, 8},
    {"_fastTest_testWilcoxonStd", (DL_FUNC) &_fastTest_testWilcoxonStd, 8},
    {"_fastTest_testWilcoxonMtx", (DL_FUNC) &_fastTest_testWilcoxonMtx, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_fastTest(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}