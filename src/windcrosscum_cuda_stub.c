/* windcrosscum_cuda_stub.c
 *
 * Fallback entry point for method="cumcuda" when the package was built
 * without CUDA support.  The real implementation is in
 * windcrosscum_cuda.cu, compiled only when HAVE_CUDA is defined
 * (see src/Makevars).
 */

#include <R.h>
#include <Rinternals.h>

#ifndef HAVE_CUDA

SEXP windcrosscum_cuda(SEXP inSeries1, SEXP inSeries2, SEXP wMax, SEXP tMax, SEXP wInc, SEXP tInc, SEXP singleS) {
    error("CUDA support was not built into this installation of wcc. "
          "Reinstall with CUDA_HOME set and nvcc available, or use method=\"cumc\".");
    return R_NilValue;
}

SEXP windcrosscum_cuda_batch(SEXP seriesArray1, SEXP seriesArray2, SEXP pairs,
                             SEXP wMax, SEXP tMax, SEXP wInc, SEXP tInc, SEXP singleS) {
    error("CUDA support was not built into this installation of wcc. "
          "Reinstall with CUDA_HOME set and nvcc available, or use method=\"cumc\".");
    return R_NilValue;
}

SEXP wccpeakpick_cuda_batch(SEXP grids, SEXP M, SEXP LsizeS, SEXP isMaxS, SEXP singleS) {
    error("CUDA support was not built into this installation of wcc. "
          "Reinstall with CUDA_HOME set and nvcc available, or use method=\"cumc\".");
    return R_NilValue;
}

SEXP wccpipeline_cuda_batch(SEXP seriesArray1, SEXP seriesArray2, SEXP pairs,
                            SEXP wMax, SEXP tMax, SEXP wInc, SEXP tInc,
                            SEXP M, SEXP LsizeS, SEXP isMaxS, SEXP singleS) {
    error("CUDA support was not built into this installation of wcc. "
          "Reinstall with CUDA_HOME set and nvcc available, or use method=\"cumc\".");
    return R_NilValue;
}

#endif
