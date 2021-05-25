// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// createModel
void createModel();
RcppExport SEXP _rfuse_createModel() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    createModel();
    return R_NilValue;
END_RCPP
}
// deleteModel
void deleteModel();
RcppExport SEXP _rfuse_deleteModel() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    deleteModel();
    return R_NilValue;
END_RCPP
}
// freeCspace
void freeCspace();
RcppExport SEXP _rfuse_freeCspace() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    freeCspace();
    return R_NilValue;
END_RCPP
}
// startTiming
void startTiming();
RcppExport SEXP _rfuse_startTiming() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    startTiming();
    return R_NilValue;
END_RCPP
}
// endTiming
void endTiming(IntegerVector print);
RcppExport SEXP _rfuse_endTiming(SEXP printSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type print(printSEXP);
    endTiming(print);
    return R_NilValue;
END_RCPP
}
// getEdgelist2
Rcpp::List getEdgelist2(Environment Sfuse, Environment Mat, IntegerVector nCORES, IntegerVector Norm2);
RcppExport SEXP _rfuse_getEdgelist2(SEXP SfuseSEXP, SEXP MatSEXP, SEXP nCORESSEXP, SEXP Norm2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type Sfuse(SfuseSEXP);
    Rcpp::traits::input_parameter< Environment >::type Mat(MatSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nCORES(nCORESSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Norm2(Norm2SEXP);
    rcpp_result_gen = Rcpp::wrap(getEdgelist2(Sfuse, Mat, nCORES, Norm2));
    return rcpp_result_gen;
END_RCPP
}
// getEdgelist
Rcpp::List getEdgelist(Environment Min, IntegerVector nCORES, IntegerVector Norm2);
RcppExport SEXP _rfuse_getEdgelist(SEXP MinSEXP, SEXP nCORESSEXP, SEXP Norm2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type Min(MinSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nCORES(nCORESSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Norm2(Norm2SEXP);
    rcpp_result_gen = Rcpp::wrap(getEdgelist(Min, nCORES, Norm2));
    return rcpp_result_gen;
END_RCPP
}
// tcpp
void tcpp(Environment Min, NumericVector nCORES);
RcppExport SEXP _rfuse_tcpp(SEXP MinSEXP, SEXP nCORESSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type Min(MinSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nCORES(nCORESSEXP);
    tcpp(Min, nCORES);
    return R_NilValue;
END_RCPP
}
// readPmatrix
void readPmatrix(Environment Pmat, NumericVector Dataset, NumericVector nCORES);
RcppExport SEXP _rfuse_readPmatrix(SEXP PmatSEXP, SEXP DatasetSEXP, SEXP nCORESSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type Pmat(PmatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Dataset(DatasetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nCORES(nCORESSEXP);
    readPmatrix(Pmat, Dataset, nCORES);
    return R_NilValue;
END_RCPP
}
// readSmatrix2
void readSmatrix2(Environment Smat, NumericVector Dataset, NumericVector nCORES);
RcppExport SEXP _rfuse_readSmatrix2(SEXP SmatSEXP, SEXP DatasetSEXP, SEXP nCORESSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type Smat(SmatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Dataset(DatasetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nCORES(nCORESSEXP);
    readSmatrix2(Smat, Dataset, nCORES);
    return R_NilValue;
END_RCPP
}
// readSmatrix
void readSmatrix(Environment Smat, NumericVector Dataset, NumericVector nCORES);
RcppExport SEXP _rfuse_readSmatrix(SEXP SmatSEXP, SEXP DatasetSEXP, SEXP nCORESSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type Smat(SmatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Dataset(DatasetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nCORES(nCORESSEXP);
    readSmatrix(Smat, Dataset, nCORES);
    return R_NilValue;
END_RCPP
}
// normalisePc
void normalisePc(Environment Pcmat, NumericVector nCORES);
RcppExport SEXP _rfuse_normalisePc(SEXP PcmatSEXP, SEXP nCORESSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type Pcmat(PcmatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nCORES(nCORESSEXP);
    normalisePc(Pcmat, nCORES);
    return R_NilValue;
END_RCPP
}
// test
void test(Environment Pmat, Environment Smat, NumericVector Dataset, NumericVector nCORES);
RcppExport SEXP _rfuse_test(SEXP PmatSEXP, SEXP SmatSEXP, SEXP DatasetSEXP, SEXP nCORESSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type Pmat(PmatSEXP);
    Rcpp::traits::input_parameter< Environment >::type Smat(SmatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Dataset(DatasetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nCORES(nCORESSEXP);
    test(Pmat, Smat, Dataset, nCORES);
    return R_NilValue;
END_RCPP
}
// findASDedges
Rcpp::List findASDedges(Environment Pfuse, IntegerVector nCORES, IntegerVector nSTEPS);
RcppExport SEXP _rfuse_findASDedges(SEXP PfuseSEXP, SEXP nCORESSEXP, SEXP nSTEPSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type Pfuse(PfuseSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nCORES(nCORESSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nSTEPS(nSTEPSSEXP);
    rcpp_result_gen = Rcpp::wrap(findASDedges(Pfuse, nCORES, nSTEPS));
    return rcpp_result_gen;
END_RCPP
}
// findCONTROLedges
Rcpp::List findCONTROLedges(Environment Pfuse, IntegerVector nCORES, IntegerVector nSTEPS);
RcppExport SEXP _rfuse_findCONTROLedges(SEXP PfuseSEXP, SEXP nCORESSEXP, SEXP nSTEPSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type Pfuse(PfuseSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nCORES(nCORESSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nSTEPS(nSTEPSSEXP);
    rcpp_result_gen = Rcpp::wrap(findCONTROLedges(Pfuse, nCORES, nSTEPS));
    return rcpp_result_gen;
END_RCPP
}
// fuseASDsets
Rcpp::List fuseASDsets(Environment Pfuse, IntegerVector nCORES, IntegerVector nTHREADS, IntegerVector nSTEPS);
RcppExport SEXP _rfuse_fuseASDsets(SEXP PfuseSEXP, SEXP nCORESSEXP, SEXP nTHREADSSEXP, SEXP nSTEPSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type Pfuse(PfuseSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nCORES(nCORESSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nTHREADS(nTHREADSSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nSTEPS(nSTEPSSEXP);
    rcpp_result_gen = Rcpp::wrap(fuseASDsets(Pfuse, nCORES, nTHREADS, nSTEPS));
    return rcpp_result_gen;
END_RCPP
}
// fuseCONTROLsets
Rcpp::List fuseCONTROLsets(Environment Pfuse, IntegerVector nCORES, IntegerVector nTHREADS, IntegerVector nSTEPS);
RcppExport SEXP _rfuse_fuseCONTROLsets(SEXP PfuseSEXP, SEXP nCORESSEXP, SEXP nTHREADSSEXP, SEXP nSTEPSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type Pfuse(PfuseSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nCORES(nCORESSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nTHREADS(nTHREADSSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nSTEPS(nSTEPSSEXP);
    rcpp_result_gen = Rcpp::wrap(fuseCONTROLsets(Pfuse, nCORES, nTHREADS, nSTEPS));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rfuse_createModel", (DL_FUNC) &_rfuse_createModel, 0},
    {"_rfuse_deleteModel", (DL_FUNC) &_rfuse_deleteModel, 0},
    {"_rfuse_freeCspace", (DL_FUNC) &_rfuse_freeCspace, 0},
    {"_rfuse_startTiming", (DL_FUNC) &_rfuse_startTiming, 0},
    {"_rfuse_endTiming", (DL_FUNC) &_rfuse_endTiming, 1},
    {"_rfuse_getEdgelist2", (DL_FUNC) &_rfuse_getEdgelist2, 4},
    {"_rfuse_getEdgelist", (DL_FUNC) &_rfuse_getEdgelist, 3},
    {"_rfuse_tcpp", (DL_FUNC) &_rfuse_tcpp, 2},
    {"_rfuse_readPmatrix", (DL_FUNC) &_rfuse_readPmatrix, 3},
    {"_rfuse_readSmatrix2", (DL_FUNC) &_rfuse_readSmatrix2, 3},
    {"_rfuse_readSmatrix", (DL_FUNC) &_rfuse_readSmatrix, 3},
    {"_rfuse_normalisePc", (DL_FUNC) &_rfuse_normalisePc, 2},
    {"_rfuse_test", (DL_FUNC) &_rfuse_test, 4},
    {"_rfuse_findASDedges", (DL_FUNC) &_rfuse_findASDedges, 3},
    {"_rfuse_findCONTROLedges", (DL_FUNC) &_rfuse_findCONTROLedges, 3},
    {"_rfuse_fuseASDsets", (DL_FUNC) &_rfuse_fuseASDsets, 4},
    {"_rfuse_fuseCONTROLsets", (DL_FUNC) &_rfuse_fuseCONTROLsets, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_rfuse(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
