#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP windcross(SEXP inSeries1, SEXP inSeries2, SEXP wMax, SEXP tMax, SEXP wInc, SEXP tInc);
extern SEXP windcrosscum(SEXP inSeries1, SEXP inSeries2, SEXP wMax, SEXP tMax, SEXP wInc, SEXP tInc);
extern SEXP windcrosscum_batch(SEXP seriesArray1, SEXP seriesArray2, SEXP pairs, SEXP wMax, SEXP tMax, SEXP wInc, SEXP tInc);
extern SEXP windcrosscum_cuda(SEXP inSeries1, SEXP inSeries2, SEXP wMax, SEXP tMax, SEXP wInc, SEXP tInc, SEXP singleS);
extern SEXP windcrosscum_cuda_batch(SEXP seriesArray1, SEXP seriesArray2, SEXP pairs, SEXP wMax, SEXP tMax, SEXP wInc, SEXP tInc, SEXP singleS);
extern SEXP wccpeakpick_cuda_batch(SEXP grids, SEXP M, SEXP LsizeS, SEXP isMaxS, SEXP singleS);
extern SEXP wccpipeline_cuda_batch(SEXP seriesArray1, SEXP seriesArray2, SEXP pairs, SEXP wMax, SEXP tMax, SEXP wInc, SEXP tInc, SEXP M, SEXP LsizeS, SEXP isMaxS, SEXP singleS);

static const R_CallMethodDef CallEntries[] = {
    {"windcross", (DL_FUNC) &windcross, 6},
    {"windcrosscum", (DL_FUNC) &windcrosscum, 6},
    {"windcrosscum_batch", (DL_FUNC) &windcrosscum_batch, 7},
    {"windcrosscum_cuda", (DL_FUNC) &windcrosscum_cuda, 7},
    {"windcrosscum_cuda_batch", (DL_FUNC) &windcrosscum_cuda_batch, 8},
    {"wccpeakpick_cuda_batch", (DL_FUNC) &wccpeakpick_cuda_batch, 5},
    {"wccpipeline_cuda_batch", (DL_FUNC) &wccpipeline_cuda_batch, 11},
    {NULL, NULL, 0}
};

void R_init_wcc(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
