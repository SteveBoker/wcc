/* windcrosscum_cuda.cu
 *
 * Batched CUDA backend for windowed cross-correlation via cumulative
 * (prefix) sums.  Same algorithm and geometry as windcrosscum.c:
 *
 *   1. per-series kernel: one block per series computes the mean, the
 *      centered series, and inclusive prefix sums of the centered values
 *      and their squares (length n+1 with leading zero) — computed once
 *      and shared by every pair that uses the series,
 *   2. fused (pair, lag-column) kernel: one block per output column scans
 *      the lagged products in shared-memory tiles (the product prefix
 *      sums never touch global memory) and evaluates all rows of the
 *      column in O(1) each from the tile scan plus the per-series prefix
 *      sums.
 *
 * All kernels are templated on the floating-point type: FP64 (default)
 * or FP32 (precision = "single" on the R side; inputs are narrowed on
 * the host, all device arithmetic runs in float, results are widened
 * back to double for R).  Missing data is not supported; the R wrapper
 * rejects NA input.
 *
 * Shared-memory tiling: the time axis is processed in tiles of TILE_N
 * samples with a running carry; the previous tile is kept in a double
 * buffer so window sums can look back up to TILE_N samples.  This
 * requires wMax <= TILE_N (checked on the host).
 *
 * Built only when the package is compiled with -DHAVE_CUDA and nvcc
 * (see src/Makevars); otherwise the stub in windcrosscum_cuda_stub.c
 * provides the entry points and raises an error.
 */

#include <cuda_runtime.h>
#include <cublas_v2.h>

#include <cstring>

#include <R.h>
#include <Rinternals.h>

#define TPB 256          /* threads per block */
#define TILE_N 2048      /* samples per shared-memory tile (TPB * 8) */
#define ITEMS_PER_THREAD (TILE_N / TPB)

#define CUDA_CHECK(call)                                                  \
    do {                                                                  \
        cudaError_t err__ = (call);                                       \
        if (err__ != cudaSuccess) {                                       \
            error("CUDA error in %s at line %d: %s", __FILE__, __LINE__,  \
                  cudaGetErrorString(err__));                             \
        }                                                                 \
    } while (0)

#define CUBLAS_CHECK(call)                                                \
    do {                                                                  \
        cublasStatus_t st__ = (call);                                     \
        if (st__ != CUBLAS_STATUS_SUCCESS) {                              \
            error("cuBLAS error in %s at line %d (status %d)", __FILE__,  \
                  __LINE__, (int) st__);                                  \
        }                                                                 \
    } while (0)

/* ------------------------------------------------------------------ */
/* Persistent per-process resources, reused across calls to amortize   */
/* cudaMalloc/cudaFree, cublasCreate, and to route transfers through   */
/* pinned memory.  All entry points run synchronously on the R main    */
/* thread, so a single set of buffers is safe.  Grow-only; released at */
/* process exit by the driver.                                         */

enum {
    BUF_IN1, BUF_IN2, BUF_XC1, BUF_XC2, BUF_CS1, BUF_CSQ1, BUF_CS2, BUF_CSQ2,
    BUF_GRID, BUF_P1, BUF_P2, BUF_M, BUF_T2, BUF_IDX, BUF_VAL, BUF_COUNT
};

static void *devBufGet(int slot, size_t bytes) {
    static void *buf[BUF_COUNT];
    static size_t cap[BUF_COUNT];
    if (cap[slot] < bytes) {
        if (buf[slot]) cudaFree(buf[slot]);
        buf[slot] = NULL;
        cap[slot] = 0;
        CUDA_CHECK(cudaMalloc(&buf[slot], bytes));
        cap[slot] = bytes;
    }
    return buf[slot];
}

/* Pinned host staging buffer for uploads/downloads. */
static void *hostStageGet(size_t bytes) {
    static void *buf = NULL;
    static size_t cap = 0;
    if (cap < bytes) {
        if (buf) cudaFreeHost(buf);
        buf = NULL;
        cap = 0;
        CUDA_CHECK(cudaHostAlloc(&buf, bytes, cudaHostAllocDefault));
        cap = bytes;
    }
    return buf;
}

static cublasHandle_t cublasHandleGet(void) {
    static cublasHandle_t handle = NULL;
    if (handle == NULL) {
        CUBLAS_CHECK(cublasCreate(&handle));
    }
    return handle;
}

template <typename T> __device__ __forceinline__ T quietNaN();
template <> __device__ __forceinline__ double quietNaN<double>() { return nan(""); }
template <> __device__ __forceinline__ float  quietNaN<float>()  { return nanf(""); }

/* Block-wide inclusive scan of vals[TILE_N] in shared memory.
 * Each thread owns ITEMS_PER_THREAD consecutive elements.  partials must
 * hold TPB elements.  After the call, vals contains the inclusive scan. */
template <typename T>
__device__ void blockScanTile(T *vals, T *partials) {
    int tid = threadIdx.x;
    int base = tid * ITEMS_PER_THREAD;
    /* serial scan of this thread's chunk */
    T sum = T(0);
    for (int i = 0; i < ITEMS_PER_THREAD; i++) {
        sum += vals[base + i];
        vals[base + i] = sum;
    }
    partials[tid] = sum;
    __syncthreads();
    /* Hillis-Steele scan of the per-thread totals */
    for (int offset = 1; offset < TPB; offset <<= 1) {
        T v = (tid >= offset) ? partials[tid - offset] : T(0);
        __syncthreads();
        if (tid >= offset) partials[tid] += v;
        __syncthreads();
    }
    /* add the exclusive prefix of preceding chunks */
    T chunkOffset = (tid > 0) ? partials[tid - 1] : T(0);
    for (int i = 0; i < ITEMS_PER_THREAD; i++) {
        vals[base + i] += chunkOffset;
    }
    __syncthreads();
}

/* One block per series.  series: input row d of a D x n column-major R
 * matrix (stride D between samples).  Outputs (contiguous per series):
 * centered values xc (length n) and prefix sums cs, csq (length n+1,
 * leading zero). */
template <typename T>
__global__ void seriesScanKernel(const T *seriesMatrix, long D, long n,
                                 T *xcAll, T *csAll, T *csqAll) {
    long d = blockIdx.x;
    int tid = threadIdx.x;
    const T *in = seriesMatrix + d;               /* sample t at in[D * t] */
    T *xc  = xcAll + d * n;
    T *cs  = csAll + d * (n + 1);
    T *csq = csqAll + d * (n + 1);

    __shared__ T partials[TPB];
    __shared__ T tile[TILE_N];
    __shared__ T carry2[2];                       /* running prefix: sum, sumsq */
    __shared__ T meanSh;

    /* pass 1: mean */
    T s = T(0);
    for (long t = tid; t < n; t += TPB) s += in[D * t];
    partials[tid] = s;
    __syncthreads();
    for (int offset = TPB >> 1; offset > 0; offset >>= 1) {
        if (tid < offset) partials[tid] += partials[tid + offset];
        __syncthreads();
    }
    if (tid == 0) {
        meanSh = partials[0] / (T) n;
        cs[0] = T(0);
        csq[0] = T(0);
        carry2[0] = T(0);
        carry2[1] = T(0);
    }
    __syncthreads();
    T mean = meanSh;

    /* pass 2: center + tiled scans of values and squares */
    for (long tileStart = 0; tileStart < n; tileStart += TILE_N) {
        /* values */
        for (int i = tid; i < TILE_N; i += TPB) {
            long t = tileStart + i;
            T v = (t < n) ? (in[D * t] - mean) : T(0);
            if (t < n) xc[t] = v;
            tile[i] = v;
        }
        __syncthreads();
        blockScanTile(tile, partials);
        for (int i = tid; i < TILE_N; i += TPB) {
            long t = tileStart + i;
            if (t < n) cs[t + 1] = tile[i] + carry2[0];
        }
        __syncthreads();
        if (tid == 0) carry2[0] += tile[TILE_N - 1];
        __syncthreads();

        /* squares */
        for (int i = tid; i < TILE_N; i += TPB) {
            long t = tileStart + i;
            T v = (t < n) ? xc[t] : T(0);
            tile[i] = v * v;
        }
        __syncthreads();
        blockScanTile(tile, partials);
        for (int i = tid; i < TILE_N; i += TPB) {
            long t = tileStart + i;
            if (t < n) csq[t + 1] = tile[i] + carry2[1];
        }
        __syncthreads();
        if (tid == 0) carry2[1] += tile[TILE_N - 1];
        __syncthreads();
    }
}

/* Fused kernel: one block per (pair, output column).
 * blockIdx.x = pair index, blockIdx.y = column index.
 * Column c maps to signed lag offset o = c - centerCol0:
 *   o >= 0: base = series1 of the pair, lagged = series2, L = o * tInc
 *   o <  0: base = series2, lagged = series1, L = -o * tInc
 * The product prefix is built tile by tile in shared memory with a
 * double buffer; window sums look back at most wMax (<= TILE_N) samples.
 * Inclusive product prefix P (P[t] = sum of products up to t): window
 * sum for prefix indices (lo, hi] is P[hi-1] - P[lo-1], with P[-1] = 0. */
template <typename T>
__global__ void wccColumnKernel(const T *xc1All, const T *cs1All, const T *csq1All,
                                const T *xc2All, const T *cs2All, const T *csq2All,
                                const int *pairs1, const int *pairs2,
                                long n, long tStart, long windowSize,
                                long windowIncrement, long lagIncrement,
                                long nRow, long nCol, long centerCol0,
                                int zeroNaN, T *out) {
    long p = blockIdx.x;
    long c = blockIdx.y;
    int tid = threadIdx.x;

    long d1 = pairs1[p] - 1;
    long d2 = pairs2[p] - 1;
    long o = c - centerCol0;
    long L = (o >= 0 ? o : -o) * lagIncrement;

    const T *base, *lagged, *csB, *csB2, *csA, *csA2;
    if (o >= 0) {
        base = xc1All + d1 * n;  csB = cs1All + d1 * (n + 1);  csB2 = csq1All + d1 * (n + 1);
        lagged = xc2All + d2 * n;  csA = cs2All + d2 * (n + 1);  csA2 = csq2All + d2 * (n + 1);
    } else {
        base = xc2All + d2 * n;  csB = cs2All + d2 * (n + 1);  csB2 = csq2All + d2 * (n + 1);
        lagged = xc1All + d1 * n;  csA = cs1All + d1 * (n + 1);  csA2 = csq1All + d1 * (n + 1);
    }
    T *outCol = out + p * nRow * nCol + c * nRow;
    T W = (T) windowSize;

    __shared__ T partials[TPB];
    __shared__ T buf[2][TILE_N];        /* current and previous tile of P */
    __shared__ T carrySh;
    if (tid == 0) carrySh = T(0);
    __syncthreads();

    int cur = 0;
    for (long tileStart = 0; tileStart < n; tileStart += TILE_N) {
        T *tile = buf[cur];
        T *prev = buf[1 - cur];
        T carry = carrySh;

        for (int i = tid; i < TILE_N; i += TPB) {
            long t = tileStart + i;
            tile[i] = (t < n && t >= L) ? base[t] * lagged[t - L] : T(0);
        }
        __syncthreads();
        blockScanTile(tile, partials);
        for (int i = tid; i < TILE_N; i += TPB) {
            tile[i] += carry;
        }
        __syncthreads();
        if (tid == 0) carrySh = tile[TILE_N - 1];
        __syncthreads();

        /* rows whose window-end prefix index hi satisfies
         * hi - 1 in [tileStart, tileStart + TILE_N - 1] */
        long kFirst = (tileStart + 1 > tStart)
                      ? (tileStart + 1 - tStart + windowIncrement - 1) / windowIncrement
                      : 0;
        long kLast = (tileStart + TILE_N - tStart) / windowIncrement;
        if (kLast > nRow - 1) kLast = nRow - 1;
        for (long k = kFirst + tid; k <= kLast; k += TPB) {
            long hi = tStart + k * windowIncrement;
            long lo = hi - windowSize;
            long hiIdx = hi - 1 - tileStart;       /* in current tile by construction */
            long loIdx = lo - 1 - tileStart;       /* may fall in previous tile */
            T Phi = tile[hiIdx];
            T Plo;
            if (lo - 1 < 0) {
                Plo = T(0);
            } else if (loIdx >= 0) {
                Plo = tile[loIdx];
            } else {
                Plo = prev[loIdx + TILE_N];
            }
            T sxy = Phi - Plo;
            T sb  = csB[hi] - csB[lo];
            T qb  = csB2[hi] - csB2[lo];
            T sa  = csA[hi - L] - csA[lo - L];
            T qa  = csA2[hi - L] - csA2[lo - L];
            T varb = W * qb - sb * sb;
            T vara = W * qa - sa * sa;
            T num = W * sxy - sb * sa;
            T den = varb * vara;
            outCol[k] = (den > T(0)) ? num / sqrt(den)
                                     : (zeroNaN ? T(0) : quietNaN<T>());
        }
        __syncthreads();
        cur = 1 - cur;
    }
}

/* ------------------------------------------------------------------ */
/* Host-side helpers: move data between R's double arrays and device   */
/* buffers of type T through the pinned staging buffer (which also     */
/* performs the FP32 narrowing/widening).                              */

static void uploadReal(double *dDst, const double *src, size_t count) {
    double *stage = (double *) hostStageGet(count * sizeof(double));
    memcpy(stage, src, count * sizeof(double));
    CUDA_CHECK(cudaMemcpy(dDst, stage, count * sizeof(double), cudaMemcpyHostToDevice));
}

static void uploadReal(float *dDst, const double *src, size_t count) {
    float *stage = (float *) hostStageGet(count * sizeof(float));
    for (size_t i = 0; i < count; i++) stage[i] = (float) src[i];
    CUDA_CHECK(cudaMemcpy(dDst, stage, count * sizeof(float), cudaMemcpyHostToDevice));
}

static void downloadReal(double *dst, const double *dSrc, size_t count) {
    double *stage = (double *) hostStageGet(count * sizeof(double));
    CUDA_CHECK(cudaMemcpy(stage, dSrc, count * sizeof(double), cudaMemcpyDeviceToHost));
    memcpy(dst, stage, count * sizeof(double));
}

static void downloadReal(double *dst, const float *dSrc, size_t count) {
    float *stage = (float *) hostStageGet(count * sizeof(float));
    CUDA_CHECK(cudaMemcpy(stage, dSrc, count * sizeof(float), cudaMemcpyDeviceToHost));
    for (size_t i = 0; i < count; i++) dst[i] = (double) stage[i];
}

/* Computes the batched WCC grid and leaves it resident on the device
 * (BUF_GRID slot, nRow x nCol column-major slabs, one per pair).
 * zeroNaN = 1 writes 0 instead of NaN for zero-variance windows (the
 * substitution the fused pipeline needs before smoothing). */
template <typename T>
static T *wccGridDevice(SEXP seriesArray1, SEXP seriesArray2, SEXP pairs,
                        SEXP wMax, SEXP tMax, SEXP wInc, SEXP tInc,
                        int zeroNaN, long *nRowOut, long *nColOut, long *nPairsOut) {
    long D1 = nrows(seriesArray1);
    long D2 = nrows(seriesArray2);
    long n = ncols(seriesArray1);
    if ((long) ncols(seriesArray2) != n) {
        error("seriesArray1 and seriesArray2 must have the same number of columns.");
    }
    long nPairs = nrows(pairs);
    const int *pairIdx = INTEGER(pairs);
    for (long p = 0; p < nPairs; p++) {
        if (pairIdx[p] < 1 || pairIdx[p] > D1 || pairIdx[p + nPairs] < 1 || pairIdx[p + nPairs] > D2) {
            error("pairs contains an index outside the series arrays.");
        }
    }

    long windowSize = (long) REAL(wMax)[0];
    long windowIncrement = (long) REAL(wInc)[0];
    long maxLag = (long) REAL(tMax)[0];
    long lagIncrement = (long) REAL(tInc)[0];

    if (windowSize > TILE_N) {
        error("method=\"cumcuda\" supports wMax up to %d; use method=\"cumc\" for larger windows.", TILE_N);
    }

    long tStart = windowSize + maxLag;
    long nRow = (n - tStart) / windowIncrement;
    long nLagSteps = maxLag / lagIncrement;
    long nCol = 2 * nLagSteps + 1;
    long centerCol0 = nLagSteps;
    if (nRow < 1) {
        error("Bad choice of parameters: the result matrix has %ld rows.", nRow);
    }
    double resCells = (double) nRow * (double) nCol * (double) nPairs;
    if (resCells > 2147483647.0) {
        error("Result array would need %.1f GB; split the pairs into smaller batches.",
              resCells * 8.0 / 1073741824.0);
    }

    T *dIn1  = (T *) devBufGet(BUF_IN1, (size_t) D1 * n * sizeof(T));
    T *dIn2  = (T *) devBufGet(BUF_IN2, (size_t) D2 * n * sizeof(T));
    T *dXc1  = (T *) devBufGet(BUF_XC1, (size_t) D1 * n * sizeof(T));
    T *dXc2  = (T *) devBufGet(BUF_XC2, (size_t) D2 * n * sizeof(T));
    T *dCs1  = (T *) devBufGet(BUF_CS1, (size_t) D1 * (n + 1) * sizeof(T));
    T *dCsq1 = (T *) devBufGet(BUF_CSQ1, (size_t) D1 * (n + 1) * sizeof(T));
    T *dCs2  = (T *) devBufGet(BUF_CS2, (size_t) D2 * (n + 1) * sizeof(T));
    T *dCsq2 = (T *) devBufGet(BUF_CSQ2, (size_t) D2 * (n + 1) * sizeof(T));
    T *dOut  = (T *) devBufGet(BUF_GRID, (size_t) nRow * nCol * nPairs * sizeof(T));
    int *dP1 = (int *) devBufGet(BUF_P1, (size_t) nPairs * sizeof(int));
    int *dP2 = (int *) devBufGet(BUF_P2, (size_t) nPairs * sizeof(int));

    uploadReal(dIn1, REAL(seriesArray1), (size_t) D1 * n);
    uploadReal(dIn2, REAL(seriesArray2), (size_t) D2 * n);
    CUDA_CHECK(cudaMemcpy(dP1, pairIdx, (size_t) nPairs * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dP2, pairIdx + nPairs, (size_t) nPairs * sizeof(int), cudaMemcpyHostToDevice));

    seriesScanKernel<<<(unsigned) D1, TPB>>>(dIn1, D1, n, dXc1, dCs1, dCsq1);
    seriesScanKernel<<<(unsigned) D2, TPB>>>(dIn2, D2, n, dXc2, dCs2, dCsq2);
    CUDA_CHECK(cudaGetLastError());

    dim3 grid((unsigned) nPairs, (unsigned) nCol);
    wccColumnKernel<<<grid, TPB>>>(dXc1, dCs1, dCsq1, dXc2, dCs2, dCsq2,
                                   dP1, dP2, n, tStart, windowSize,
                                   windowIncrement, lagIncrement,
                                   nRow, nCol, centerCol0, zeroNaN, dOut);
    CUDA_CHECK(cudaGetLastError());

    *nRowOut = nRow;
    *nColOut = nCol;
    *nPairsOut = nPairs;
    return dOut;
}

template <typename T>
static SEXP windcrosscum_cuda_batch_impl(SEXP seriesArray1, SEXP seriesArray2, SEXP pairs,
                                         SEXP wMax, SEXP tMax, SEXP wInc, SEXP tInc) {
    long nRow, nCol, nPairs;
    T *dOut = wccGridDevice<T>(seriesArray1, seriesArray2, pairs,
                               wMax, tMax, wInc, tInc, 0,
                               &nRow, &nCol, &nPairs);

    SEXP dims = PROTECT(allocVector(INTSXP, 3));
    INTEGER(dims)[0] = (int) nRow;
    INTEGER(dims)[1] = (int) nCol;
    INTEGER(dims)[2] = (int) nPairs;
    SEXP corResult = PROTECT(allocArray(REALSXP, dims));
    double *out = REAL(corResult);

    downloadReal(out, dOut, (size_t) nRow * nCol * nPairs);

    /* NaN (zero-variance windows) -> R NA */
    long total = nRow * nCol * nPairs;
    for (long j = 0; j < total; j++) {
        if (isnan(out[j])) out[j] = NA_REAL;
    }

    UNPROTECT(2);
    return corResult;
}

static int precisionIsSingle(SEXP singleS) {
    if (!isInteger(singleS) && !isLogical(singleS) && !isReal(singleS)) {
        error("precision flag must be a scalar.");
    }
    return asInteger(singleS) != 0;
}

extern "C" SEXP windcrosscum_cuda_batch(SEXP seriesArray1, SEXP seriesArray2, SEXP pairs,
                                        SEXP wMax, SEXP tMax, SEXP wInc, SEXP tInc,
                                        SEXP singleS) {
    if (!isReal(seriesArray1) || !isReal(seriesArray2) || !isMatrix(seriesArray1) || !isMatrix(seriesArray2)) {
        error("seriesArray1 and seriesArray2 must be numeric matrices.");
    }
    if (!isInteger(pairs) || !isMatrix(pairs) || ncols(pairs) != 2) {
        error("pairs must be an integer matrix with two columns.");
    }
    if (precisionIsSingle(singleS)) {
        return windcrosscum_cuda_batch_impl<float>(seriesArray1, seriesArray2, pairs,
                                                   wMax, tMax, wInc, tInc);
    }
    return windcrosscum_cuda_batch_impl<double>(seriesArray1, seriesArray2, pairs,
                                                wMax, tMax, wInc, tInc);
}

/* Single-dyad entry point: shim over the batch routine with one (1, 1) pair. */
extern "C" SEXP windcrosscum_cuda(SEXP inSeries1, SEXP inSeries2, SEXP wMax,
                                  SEXP tMax, SEXP wInc, SEXP tInc, SEXP singleS) {
    if (!isReal(inSeries1) || !isReal(inSeries2)) {
        error("inSeries1 and inSeries2 must be numeric.");
    }
    R_xlen_t n = XLENGTH(inSeries1);
    if (XLENGTH(inSeries2) != n) {
        error("inSeries1 and inSeries2 must have equal length.");
    }

    SEXP m1 = PROTECT(allocMatrix(REALSXP, 1, (int) n));
    SEXP m2 = PROTECT(allocMatrix(REALSXP, 1, (int) n));
    memcpy(REAL(m1), REAL(inSeries1), n * sizeof(double));
    memcpy(REAL(m2), REAL(inSeries2), n * sizeof(double));
    SEXP onePair = PROTECT(allocMatrix(INTSXP, 1, 2));
    INTEGER(onePair)[0] = 1;
    INTEGER(onePair)[1] = 1;

    SEXP res3d = PROTECT(windcrosscum_cuda_batch(m1, m2, onePair, wMax, tMax, wInc, tInc, singleS));

    int *d = INTEGER(getAttrib(res3d, R_DimSymbol));
    SEXP corResult = PROTECT(allocMatrix(REALSXP, d[0], d[1]));
    memcpy(REAL(corResult), REAL(res3d), (size_t) d[0] * d[1] * sizeof(double));

    UNPROTECT(5);
    return corResult;
}

/* ------------------------------------------------------------------ */
/* Batched peak picking.                                               */
/*                                                                     */
/* The loess+spline smoother is applied as a precomputed linear        */
/* operator M ((2*colLen-1) x colLen, built and cached on the R side   */
/* by wccSmoothMatrix): all P grids are smoothed with one strided-     */
/* batched GEMM, T2[p] = grid[p] %*% t(M), then one thread per         */
/* (grid, row) runs the expanding-window look-ahead search and writes  */
/* peak index and value.  Only the P x nRow index/value arrays are     */
/* downloaded.                                                         */

/* One thread per (grid, row).  t2 row r of slab p is read column-wise
 * from the GEMM output (column-major nRow x m per slab).  Replicates
 * the R search: mx[j] = max over half-width-j window around the center
 * (0-based c0 = colLen-1); strict improvements reset the look-ahead,
 * Lsize consecutive non-improvements stop the search; first index in
 * the final window matching the running max gives the peak position;
 * positions beyond colLen - Lsize - 1 fail to NaN. */
template <typename T>
__global__ void peakSearchKernel(const T *T2, long nRow, long colLen,
                                 long P, long Lsize, int findMax,
                                 T *outIndex, T *outValue) {
    long g = (long) blockIdx.x * blockDim.x + threadIdx.x;
    if (g >= nRow * P) return;
    long p = g / nRow;
    long r = g % nRow;
    long m = 2 * colLen - 1;
    const T *slab = T2 + p * nRow * m;
    long c0 = colLen - 1;
    T sign = findMax ? T(1) : T(-1);

#define T2AT(j) (sign * slab[r + nRow * (j)])

    /* expanding-window maxima with look-ahead, computed on the fly:
     * window half-width j keeps running maxima of the left and right
     * arms, so mx[j] needs only two new reads per step. */
    T lv = T2AT(c0);
    T rv = lv;
    T mmx = T(0);
    long lookAhead = 0;
    long windowWidth = colLen - 1;
    for (long j = 1; j <= colLen - 1; j++) {
        T a = T2AT(c0 - j);
        T b = T2AT(c0 + j);
        if (a > lv) lv = a;
        if (b > rv) rv = b;
        T mx = (lv > rv) ? lv : rv;
        if (j == 1) {
            mmx = mx;
        } else if (mx > mmx) {
            lookAhead = 0;
            mmx = mx;
        } else {
            lookAhead++;
            if (lookAhead >= Lsize) {
                windowWidth = j;
                break;
            }
        }
    }

    /* first match of mmx within the final window */
    long index = -1;
    for (long t = c0 - windowWidth; t <= c0 + windowWidth; t++) {
        if (T2AT(t) == mmx) {
            index = t;
            break;
        }
    }
    long position = index - c0;
    if (position > (colLen - Lsize - 1) || position < -(colLen - Lsize - 1)) {
        outIndex[g] = quietNaN<T>();
        outValue[g] = quietNaN<T>();
    } else {
        outIndex[g] = (T) position;
        outValue[g] = sign * mmx;
    }
#undef T2AT
}

static cublasStatus_t gemmStridedBatched(cublasHandle_t handle,
                                         int nRow, int m, int colLen,
                                         const double *dGrids, const double *dM,
                                         double *dT2, int P) {
    const double one = 1.0, zero = 0.0;
    return cublasDgemmStridedBatched(handle,
        CUBLAS_OP_N, CUBLAS_OP_T,
        nRow, m, colLen,
        &one,
        dGrids, nRow, (long long) nRow * colLen,
        dM, m, 0,
        &zero,
        dT2, nRow, (long long) nRow * m,
        P);
}

static cublasStatus_t gemmStridedBatched(cublasHandle_t handle,
                                         int nRow, int m, int colLen,
                                         const float *dGrids, const float *dM,
                                         float *dT2, int P) {
    const float one = 1.0f, zero = 0.0f;
    return cublasSgemmStridedBatched(handle,
        CUBLAS_OP_N, CUBLAS_OP_T,
        nRow, m, colLen,
        &one,
        dGrids, nRow, (long long) nRow * colLen,
        dM, m, 0,
        &zero,
        dT2, nRow, (long long) nRow * m,
        P);
}

/* Smooths + peak-searches P device-resident grids (dGrids, nRow x
 * colLen column-major slabs).  Uploads M, runs the batched GEMM and the
 * search kernel, downloads only the P x nRow index/value arrays, and
 * returns list(index, value) with NaN -> NA. */
template <typename T>
static SEXP peakPickDevice(T *dGrids, SEXP M, long nRow, long colLen,
                           long P, long Lsize, int findMax) {
    long m = 2 * colLen - 1;
    double t2Cells = (double) nRow * (double) m * (double) P;
    if (t2Cells > 2147483647.0) {
        error("Smoothed array would need %.1f GB; split the pairs into smaller batches.",
              t2Cells * 8.0 / 1073741824.0);
    }

    SEXP outIndexS = PROTECT(allocVector(REALSXP, nRow * P));
    SEXP outValueS = PROTECT(allocVector(REALSXP, nRow * P));

    T *dM   = (T *) devBufGet(BUF_M, (size_t) m * colLen * sizeof(T));
    T *dT2  = (T *) devBufGet(BUF_T2, (size_t) nRow * m * P * sizeof(T));
    T *dIdx = (T *) devBufGet(BUF_IDX, (size_t) nRow * P * sizeof(T));
    T *dVal = (T *) devBufGet(BUF_VAL, (size_t) nRow * P * sizeof(T));

    uploadReal(dM, REAL(M), (size_t) m * colLen);

    /* T2[p] (nRow x m) = grid[p] (nRow x colLen) * M^T (colLen x m) */
    CUBLAS_CHECK(gemmStridedBatched(cublasHandleGet(), (int) nRow, (int) m, (int) colLen,
                                    dGrids, dM, dT2, (int) P));

    long total = nRow * P;
    long nBlocks = (total + TPB - 1) / TPB;
    peakSearchKernel<<<(unsigned) nBlocks, TPB>>>(dT2, nRow, colLen, P, Lsize, findMax,
                                                  dIdx, dVal);
    CUDA_CHECK(cudaGetLastError());

    downloadReal(REAL(outIndexS), dIdx, (size_t) total);
    downloadReal(REAL(outValueS), dVal, (size_t) total);

    double *oi = REAL(outIndexS);
    double *ov = REAL(outValueS);
    for (long j = 0; j < total; j++) {
        if (isnan(oi[j])) { oi[j] = NA_REAL; ov[j] = NA_REAL; }
    }

    SEXP res = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(res, 0, outIndexS);
    SET_VECTOR_ELT(res, 1, outValueS);
    UNPROTECT(3);
    return res;
}

template <typename T>
static SEXP wccpeakpick_cuda_batch_impl(SEXP grids, SEXP M, long nRow, long colLen,
                                        long P, long Lsize, int findMax) {
    T *dGrids = (T *) devBufGet(BUF_GRID, (size_t) nRow * colLen * P * sizeof(T));
    uploadReal(dGrids, REAL(grids), (size_t) nRow * colLen * P);
    return peakPickDevice<T>(dGrids, M, nRow, colLen, P, Lsize, findMax);
}

/* grids: nRow x colLen x P array; M: (2*colLen-1) x colLen smoother
 * matrix; LsizeS, isMaxS: scalar numerics; singleS: 1 = FP32, 0 = FP64.
 * Returns list(index, value), each a numeric vector of length
 * nRow * P (NaN -> NA). */
extern "C" SEXP wccpeakpick_cuda_batch(SEXP grids, SEXP M, SEXP LsizeS, SEXP isMaxS,
                                       SEXP singleS) {
    if (!isReal(grids) || !isReal(M)) {
        error("grids and M must be numeric.");
    }
    SEXP gDims = getAttrib(grids, R_DimSymbol);
    if (isNull(gDims) || LENGTH(gDims) != 3) {
        error("grids must be a 3-dimensional array.");
    }
    long nRow = INTEGER(gDims)[0];
    long colLen = INTEGER(gDims)[1];
    long P = INTEGER(gDims)[2];
    long m = 2 * colLen - 1;
    if (nrows(M) != m || ncols(M) != colLen) {
        error("M must be a (2*colLen-1) x colLen matrix.");
    }
    long Lsize = (long) REAL(LsizeS)[0];
    int findMax = (int) REAL(isMaxS)[0];

    if (precisionIsSingle(singleS)) {
        return wccpeakpick_cuda_batch_impl<float>(grids, M, nRow, colLen, P, Lsize, findMax);
    }
    return wccpeakpick_cuda_batch_impl<double>(grids, M, nRow, colLen, P, Lsize, findMax);
}

/* ------------------------------------------------------------------ */
/* Fused pipeline: WCC + peak pick without the grid ever leaving the   */
/* device.  Zero-variance cells are written as 0 (matching the         */
/* g[is.na(g)] <- 0 substitution of the two-step pipeline), the grid   */
/* stays resident in BUF_GRID, and only the P x nRow index/value       */
/* arrays cross back to the host.                                      */

template <typename T>
static SEXP wccpipeline_cuda_batch_impl(SEXP seriesArray1, SEXP seriesArray2, SEXP pairs,
                                        SEXP wMax, SEXP tMax, SEXP wInc, SEXP tInc,
                                        SEXP M, long Lsize, int findMax) {
    long nRow, nCol, nPairs;
    T *dGrid = wccGridDevice<T>(seriesArray1, seriesArray2, pairs,
                                wMax, tMax, wInc, tInc, 1,
                                &nRow, &nCol, &nPairs);
    if (nrows(M) != 2 * nCol - 1 || ncols(M) != nCol) {
        error("M must be a (2*nCol-1) x nCol matrix.");
    }
    return peakPickDevice<T>(dGrid, M, nRow, nCol, nPairs, Lsize, findMax);
}

/* seriesArray1/2: D x n double matrices; pairs: P x 2 integer matrix;
 * M: (2*nCol-1) x nCol smoother matrix from wccSmoothMatrix; LsizeS,
 * isMaxS: scalar numerics; singleS: 1 = FP32, 0 = FP64.
 * Returns list(index, value), each of length nRow * P (NaN -> NA). */
extern "C" SEXP wccpipeline_cuda_batch(SEXP seriesArray1, SEXP seriesArray2, SEXP pairs,
                                       SEXP wMax, SEXP tMax, SEXP wInc, SEXP tInc,
                                       SEXP M, SEXP LsizeS, SEXP isMaxS, SEXP singleS) {
    if (!isReal(seriesArray1) || !isReal(seriesArray2) || !isMatrix(seriesArray1) || !isMatrix(seriesArray2)) {
        error("seriesArray1 and seriesArray2 must be numeric matrices.");
    }
    if (!isInteger(pairs) || !isMatrix(pairs) || ncols(pairs) != 2) {
        error("pairs must be an integer matrix with two columns.");
    }
    if (!isReal(M) || !isMatrix(M)) {
        error("M must be a numeric matrix.");
    }
    long Lsize = (long) REAL(LsizeS)[0];
    int findMax = (int) REAL(isMaxS)[0];

    if (precisionIsSingle(singleS)) {
        return wccpipeline_cuda_batch_impl<float>(seriesArray1, seriesArray2, pairs,
                                                  wMax, tMax, wInc, tInc, M, Lsize, findMax);
    }
    return wccpipeline_cuda_batch_impl<double>(seriesArray1, seriesArray2, pairs,
                                               wMax, tMax, wInc, tInc, M, Lsize, findMax);
}
