// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// chopthin
List chopthin(std::vector<double>& w, int N, double eta);
RcppExport SEXP chopthin_chopthin(SEXP wSEXP, SEXP NSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::vector<double>& >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    __result = Rcpp::wrap(chopthin(w, N, eta));
    return __result;
END_RCPP
}
