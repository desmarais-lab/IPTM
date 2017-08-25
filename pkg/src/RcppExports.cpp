// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// callRMultinom
IntegerVector callRMultinom(NumericVector x);
RcppExport SEXP _IPTM_callRMultinom(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(callRMultinom(x));
    return rcpp_result_gen;
END_RCPP
}
// multinom_vec
IntegerVector multinom_vec(int nSample, NumericVector props);
RcppExport SEXP _IPTM_multinom_vec(SEXP nSampleSEXP, SEXP propsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nSample(nSampleSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type props(propsSEXP);
    rcpp_result_gen = Rcpp::wrap(multinom_vec(nSample, props));
    return rcpp_result_gen;
END_RCPP
}
// which_int
int which_int(int value, IntegerVector x);
RcppExport SEXP _IPTM_which_int(SEXP valueSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type value(valueSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(which_int(value, x));
    return rcpp_result_gen;
END_RCPP
}
// which_num
int which_num(int value, NumericVector x);
RcppExport SEXP _IPTM_which_num(SEXP valueSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type value(valueSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(which_num(value, x));
    return rcpp_result_gen;
END_RCPP
}
// rdirichlet_cpp
arma::mat rdirichlet_cpp(int num_samples, arma::vec alpha_m);
RcppExport SEXP _IPTM_rdirichlet_cpp(SEXP num_samplesSEXP, SEXP alpha_mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type num_samples(num_samplesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha_m(alpha_mSEXP);
    rcpp_result_gen = Rcpp::wrap(rdirichlet_cpp(num_samples, alpha_m));
    return rcpp_result_gen;
END_RCPP
}
// rbinom_mat
IntegerMatrix rbinom_mat(NumericMatrix probmat);
RcppExport SEXP _IPTM_rbinom_mat(SEXP probmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type probmat(probmatSEXP);
    rcpp_result_gen = Rcpp::wrap(rbinom_mat(probmat));
    return rcpp_result_gen;
END_RCPP
}
// which_cpp
IntegerVector which_cpp(int value, NumericVector x);
RcppExport SEXP _IPTM_which_cpp(SEXP valueSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type value(valueSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(which_cpp(value, x));
    return rcpp_result_gen;
END_RCPP
}
// pdmat
NumericMatrix pdmat(List currentZ, NumericVector currentC, int nIP);
RcppExport SEXP _IPTM_pdmat(SEXP currentZSEXP, SEXP currentCSEXP, SEXP nIPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type currentZ(currentZSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type currentC(currentCSEXP);
    Rcpp::traits::input_parameter< int >::type nIP(nIPSEXP);
    rcpp_result_gen = Rcpp::wrap(pdmat(currentZ, currentC, nIP));
    return rcpp_result_gen;
END_RCPP
}
// History
List History(List edge, NumericMatrix p_d, IntegerVector node, double when);
RcppExport SEXP _IPTM_History(SEXP edgeSEXP, SEXP p_dSEXP, SEXP nodeSEXP, SEXP whenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type edge(edgeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type p_d(p_dSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type node(nodeSEXP);
    Rcpp::traits::input_parameter< double >::type when(whenSEXP);
    rcpp_result_gen = Rcpp::wrap(History(edge, p_d, node, when));
    return rcpp_result_gen;
END_RCPP
}
// Degree
List Degree(List history, IntegerVector node, int sender);
RcppExport SEXP _IPTM_Degree(SEXP historySEXP, SEXP nodeSEXP, SEXP senderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type history(historySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type node(nodeSEXP);
    Rcpp::traits::input_parameter< int >::type sender(senderSEXP);
    rcpp_result_gen = Rcpp::wrap(Degree(history, node, sender));
    return rcpp_result_gen;
END_RCPP
}
// Dyadic
List Dyadic(List history, IntegerVector node, int sender);
RcppExport SEXP _IPTM_Dyadic(SEXP historySEXP, SEXP nodeSEXP, SEXP senderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type history(historySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type node(nodeSEXP);
    Rcpp::traits::input_parameter< int >::type sender(senderSEXP);
    rcpp_result_gen = Rcpp::wrap(Dyadic(history, node, sender));
    return rcpp_result_gen;
END_RCPP
}
// Triadic
List Triadic(List history, IntegerVector node, int sender);
RcppExport SEXP _IPTM_Triadic(SEXP historySEXP, SEXP nodeSEXP, SEXP senderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type history(historySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type node(nodeSEXP);
    Rcpp::traits::input_parameter< int >::type sender(senderSEXP);
    rcpp_result_gen = Rcpp::wrap(Triadic(history, node, sender));
    return rcpp_result_gen;
END_RCPP
}
// Triadic_reduced
List Triadic_reduced(List triadic);
RcppExport SEXP _IPTM_Triadic_reduced(SEXP triadicSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type triadic(triadicSEXP);
    rcpp_result_gen = Rcpp::wrap(Triadic_reduced(triadic));
    return rcpp_result_gen;
END_RCPP
}
// Netstats_cpp
List Netstats_cpp(List historyIP, IntegerVector node, IntegerVector netstat);
RcppExport SEXP _IPTM_Netstats_cpp(SEXP historyIPSEXP, SEXP nodeSEXP, SEXP netstatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type historyIP(historyIPSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type node(nodeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type netstat(netstatSEXP);
    rcpp_result_gen = Rcpp::wrap(Netstats_cpp(historyIP, node, netstat));
    return rcpp_result_gen;
END_RCPP
}
// MultiplyXB
NumericVector MultiplyXB(NumericMatrix X, NumericVector B);
RcppExport SEXP _IPTM_MultiplyXB(SEXP XSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(MultiplyXB(X, B));
    return rcpp_result_gen;
END_RCPP
}
// MultiplyXBList
List MultiplyXBList(List X, List B);
RcppExport SEXP _IPTM_MultiplyXBList(SEXP XSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type X(XSEXP);
    Rcpp::traits::input_parameter< List >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(MultiplyXBList(X, B));
    return rcpp_result_gen;
END_RCPP
}
// UpdateDenom
double UpdateDenom(double alpha, IntegerVector nwordtable);
RcppExport SEXP _IPTM_UpdateDenom(SEXP alphaSEXP, SEXP nwordtableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nwordtable(nwordtableSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateDenom(alpha, nwordtable));
    return rcpp_result_gen;
END_RCPP
}
// UpdateNum
NumericVector UpdateNum(NumericVector vec, List nKwordtable);
RcppExport SEXP _IPTM_UpdateNum(SEXP vecSEXP, SEXP nKwordtableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< List >::type nKwordtable(nKwordtableSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateNum(vec, nKwordtable));
    return rcpp_result_gen;
END_RCPP
}
// tabulateC
IntegerVector tabulateC(const IntegerVector& x, const signed max);
RcppExport SEXP _IPTM_tabulateC(SEXP xSEXP, SEXP maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const signed >::type max(maxSEXP);
    rcpp_result_gen = Rcpp::wrap(tabulateC(x, max));
    return rcpp_result_gen;
END_RCPP
}
// lambda_cpp
arma::mat lambda_cpp(arma::vec p_d, List XB);
RcppExport SEXP _IPTM_lambda_cpp(SEXP p_dSEXP, SEXP XBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type p_d(p_dSEXP);
    Rcpp::traits::input_parameter< List >::type XB(XBSEXP);
    rcpp_result_gen = Rcpp::wrap(lambda_cpp(p_d, XB));
    return rcpp_result_gen;
END_RCPP
}
// TopicInEqZ
NumericVector TopicInEqZ(int K, IntegerVector currentZ_d, double alpha, NumericVector mvec);
RcppExport SEXP _IPTM_TopicInEqZ(SEXP KSEXP, SEXP currentZ_dSEXP, SEXP alphaSEXP, SEXP mvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type currentZ_d(currentZ_dSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mvec(mvecSEXP);
    rcpp_result_gen = Rcpp::wrap(TopicInEqZ(K, currentZ_d, alpha, mvec));
    return rcpp_result_gen;
END_RCPP
}
// WordInEqZ
NumericMatrix WordInEqZ(int K, IntegerVector textlistd, List tableW, double beta, NumericVector nvec);
RcppExport SEXP _IPTM_WordInEqZ(SEXP KSEXP, SEXP textlistdSEXP, SEXP tableWSEXP, SEXP betaSEXP, SEXP nvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type textlistd(textlistdSEXP);
    Rcpp::traits::input_parameter< List >::type tableW(tableWSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nvec(nvecSEXP);
    rcpp_result_gen = Rcpp::wrap(WordInEqZ(K, textlistd, tableW, beta, nvec));
    return rcpp_result_gen;
END_RCPP
}
// EdgeInEqZ
double EdgeInEqZ(IntegerMatrix iJi, NumericMatrix lambda, double delta);
RcppExport SEXP _IPTM_EdgeInEqZ(SEXP iJiSEXP, SEXP lambdaSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type iJi(iJiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(EdgeInEqZ(iJi, lambda, delta));
    return rcpp_result_gen;
END_RCPP
}
// EdgeInEqZ_Gibbs
double EdgeInEqZ_Gibbs(arma::mat iJi, arma::mat lambda, double delta);
RcppExport SEXP _IPTM_EdgeInEqZ_Gibbs(SEXP iJiSEXP, SEXP lambdaSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type iJi(iJiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(EdgeInEqZ_Gibbs(iJi, lambda, delta));
    return rcpp_result_gen;
END_RCPP
}
// EdgeInEqZ_Gibbs2
arma::vec EdgeInEqZ_Gibbs2(arma::mat iJi, arma::mat lambda, double delta);
RcppExport SEXP _IPTM_EdgeInEqZ_Gibbs2(SEXP iJiSEXP, SEXP lambdaSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type iJi(iJiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(EdgeInEqZ_Gibbs2(iJi, lambda, delta));
    return rcpp_result_gen;
END_RCPP
}
// TimeObsInEqZ
double TimeObsInEqZ(NumericVector LambdaiJi, double observedtdiff, double observediJi);
RcppExport SEXP _IPTM_TimeObsInEqZ(SEXP LambdaiJiSEXP, SEXP observedtdiffSEXP, SEXP observediJiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type LambdaiJi(LambdaiJiSEXP);
    Rcpp::traits::input_parameter< double >::type observedtdiff(observedtdiffSEXP);
    Rcpp::traits::input_parameter< double >::type observediJi(observediJiSEXP);
    rcpp_result_gen = Rcpp::wrap(TimeObsInEqZ(LambdaiJi, observedtdiff, observediJi));
    return rcpp_result_gen;
END_RCPP
}
// lambdaiJi
NumericVector lambdaiJi(NumericVector p_d, List XB, IntegerMatrix iJi);
RcppExport SEXP _IPTM_lambdaiJi(SEXP p_dSEXP, SEXP XBSEXP, SEXP iJiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p_d(p_dSEXP);
    Rcpp::traits::input_parameter< List >::type XB(XBSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type iJi(iJiSEXP);
    rcpp_result_gen = Rcpp::wrap(lambdaiJi(p_d, XB, iJi));
    return rcpp_result_gen;
END_RCPP
}
// DataAug_cpp
arma::vec DataAug_cpp(arma::vec iJi_di, arma::vec lambda_di, List XB, arma::vec p_d, double delta, double timeinc_d, int i, int j);
RcppExport SEXP _IPTM_DataAug_cpp(SEXP iJi_diSEXP, SEXP lambda_diSEXP, SEXP XBSEXP, SEXP p_dSEXP, SEXP deltaSEXP, SEXP timeinc_dSEXP, SEXP iSEXP, SEXP jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type iJi_di(iJi_diSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda_di(lambda_diSEXP);
    Rcpp::traits::input_parameter< List >::type XB(XBSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p_d(p_dSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< double >::type timeinc_d(timeinc_dSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    rcpp_result_gen = Rcpp::wrap(DataAug_cpp(iJi_di, lambda_di, XB, p_d, delta, timeinc_d, i, j));
    return rcpp_result_gen;
END_RCPP
}
// DataAug_cpp_Gibbs
arma::vec DataAug_cpp_Gibbs(arma::vec iJi_di, arma::vec lambda_di, List XB, arma::vec p_d, double delta, double timeinc_d, int j);
RcppExport SEXP _IPTM_DataAug_cpp_Gibbs(SEXP iJi_diSEXP, SEXP lambda_diSEXP, SEXP XBSEXP, SEXP p_dSEXP, SEXP deltaSEXP, SEXP timeinc_dSEXP, SEXP jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type iJi_di(iJi_diSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda_di(lambda_diSEXP);
    Rcpp::traits::input_parameter< List >::type XB(XBSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p_d(p_dSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< double >::type timeinc_d(timeinc_dSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    rcpp_result_gen = Rcpp::wrap(DataAug_cpp_Gibbs(iJi_di, lambda_di, XB, p_d, delta, timeinc_d, j));
    return rcpp_result_gen;
END_RCPP
}
// DataAug_cpp_Gibbs_noObs
arma::vec DataAug_cpp_Gibbs_noObs(arma::vec iJi_di, arma::vec lambda_di, List XB, arma::vec p_d, double delta, int j);
RcppExport SEXP _IPTM_DataAug_cpp_Gibbs_noObs(SEXP iJi_diSEXP, SEXP lambda_diSEXP, SEXP XBSEXP, SEXP p_dSEXP, SEXP deltaSEXP, SEXP jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type iJi_di(iJi_diSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda_di(lambda_diSEXP);
    Rcpp::traits::input_parameter< List >::type XB(XBSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p_d(p_dSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    rcpp_result_gen = Rcpp::wrap(DataAug_cpp_Gibbs_noObs(iJi_di, lambda_di, XB, p_d, delta, j));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_IPTM();

static const R_CallMethodDef CallEntries[] = {
    {"_IPTM_callRMultinom", (DL_FUNC) &_IPTM_callRMultinom, 1},
    {"_IPTM_multinom_vec", (DL_FUNC) &_IPTM_multinom_vec, 2},
    {"_IPTM_which_int", (DL_FUNC) &_IPTM_which_int, 2},
    {"_IPTM_which_num", (DL_FUNC) &_IPTM_which_num, 2},
    {"_IPTM_rdirichlet_cpp", (DL_FUNC) &_IPTM_rdirichlet_cpp, 2},
    {"_IPTM_rbinom_mat", (DL_FUNC) &_IPTM_rbinom_mat, 1},
    {"_IPTM_which_cpp", (DL_FUNC) &_IPTM_which_cpp, 2},
    {"_IPTM_pdmat", (DL_FUNC) &_IPTM_pdmat, 3},
    {"_IPTM_History", (DL_FUNC) &_IPTM_History, 4},
    {"_IPTM_Degree", (DL_FUNC) &_IPTM_Degree, 3},
    {"_IPTM_Dyadic", (DL_FUNC) &_IPTM_Dyadic, 3},
    {"_IPTM_Triadic", (DL_FUNC) &_IPTM_Triadic, 3},
    {"_IPTM_Triadic_reduced", (DL_FUNC) &_IPTM_Triadic_reduced, 1},
    {"_IPTM_Netstats_cpp", (DL_FUNC) &_IPTM_Netstats_cpp, 3},
    {"_IPTM_MultiplyXB", (DL_FUNC) &_IPTM_MultiplyXB, 2},
    {"_IPTM_MultiplyXBList", (DL_FUNC) &_IPTM_MultiplyXBList, 2},
    {"_IPTM_UpdateDenom", (DL_FUNC) &_IPTM_UpdateDenom, 2},
    {"_IPTM_UpdateNum", (DL_FUNC) &_IPTM_UpdateNum, 2},
    {"_IPTM_tabulateC", (DL_FUNC) &_IPTM_tabulateC, 2},
    {"_IPTM_lambda_cpp", (DL_FUNC) &_IPTM_lambda_cpp, 2},
    {"_IPTM_TopicInEqZ", (DL_FUNC) &_IPTM_TopicInEqZ, 4},
    {"_IPTM_WordInEqZ", (DL_FUNC) &_IPTM_WordInEqZ, 5},
    {"_IPTM_EdgeInEqZ", (DL_FUNC) &_IPTM_EdgeInEqZ, 3},
    {"_IPTM_EdgeInEqZ_Gibbs", (DL_FUNC) &_IPTM_EdgeInEqZ_Gibbs, 3},
    {"_IPTM_EdgeInEqZ_Gibbs2", (DL_FUNC) &_IPTM_EdgeInEqZ_Gibbs2, 3},
    {"_IPTM_TimeObsInEqZ", (DL_FUNC) &_IPTM_TimeObsInEqZ, 3},
    {"_IPTM_lambdaiJi", (DL_FUNC) &_IPTM_lambdaiJi, 3},
    {"_IPTM_DataAug_cpp", (DL_FUNC) &_IPTM_DataAug_cpp, 8},
    {"_IPTM_DataAug_cpp_Gibbs", (DL_FUNC) &_IPTM_DataAug_cpp_Gibbs, 7},
    {"_IPTM_DataAug_cpp_Gibbs_noObs", (DL_FUNC) &_IPTM_DataAug_cpp_Gibbs_noObs, 6},
    {"_rcpp_module_boot_IPTM", (DL_FUNC) &_rcpp_module_boot_IPTM, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_IPTM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
