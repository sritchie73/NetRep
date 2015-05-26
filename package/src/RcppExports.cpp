// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// AdjProps
List AdjProps(SEXP pAdjacency, IntegerVector subsetIndices);
RcppExport SEXP fastModPres_AdjProps(SEXP pAdjacencySEXP, SEXP subsetIndicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type pAdjacency(pAdjacencySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type subsetIndices(subsetIndicesSEXP);
    __result = Rcpp::wrap(AdjProps(pAdjacency, subsetIndices));
    return __result;
END_RCPP
}
// CoexpStats
List CoexpStats(SEXP pCoexpD, IntegerVector discIndices, SEXP pCoexpT, IntegerVector testIndices);
RcppExport SEXP fastModPres_CoexpStats(SEXP pCoexpDSEXP, SEXP discIndicesSEXP, SEXP pCoexpTSEXP, SEXP testIndicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type pCoexpD(pCoexpDSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type discIndices(discIndicesSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pCoexpT(pCoexpTSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type testIndices(testIndicesSEXP);
    __result = Rcpp::wrap(CoexpStats(pCoexpD, discIndices, pCoexpT, testIndices));
    return __result;
END_RCPP
}
// RangeSubset
List RangeSubset(SEXP pDat, IntegerVector subsetIndices);
RcppExport SEXP fastModPres_RangeSubset(SEXP pDatSEXP, SEXP subsetIndicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type pDat(pDatSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type subsetIndices(subsetIndicesSEXP);
    __result = Rcpp::wrap(RangeSubset(pDat, subsetIndices));
    return __result;
END_RCPP
}
// BigRange
List BigRange(SEXP pDat);
RcppExport SEXP fastModPres_BigRange(SEXP pDatSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type pDat(pDatSEXP);
    __result = Rcpp::wrap(BigRange(pDat));
    return __result;
END_RCPP
}
// CheckFinite
void CheckFinite(SEXP pDat);
RcppExport SEXP fastModPres_CheckFinite(SEXP pDatSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type pDat(pDatSEXP);
    CheckFinite(pDat);
    return R_NilValue;
END_RCPP
}
// DataProps
List DataProps(SEXP pDat, IntegerVector subsetIndices);
RcppExport SEXP fastModPres_DataProps(SEXP pDatSEXP, SEXP subsetIndicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type pDat(pDatSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type subsetIndices(subsetIndicesSEXP);
    __result = Rcpp::wrap(DataProps(pDat, subsetIndices));
    return __result;
END_RCPP
}
// Scale
void Scale(SEXP pDat, SEXP spDat);
RcppExport SEXP fastModPres_Scale(SEXP pDatSEXP, SEXP spDatSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type pDat(pDatSEXP);
    Rcpp::traits::input_parameter< SEXP >::type spDat(spDatSEXP);
    Scale(pDat, spDat);
    return R_NilValue;
END_RCPP
}
