# ---------------------------------------------------------------------
# Program: testWccCalcMethods.R
#
# Correctness tests for the wccCalc backends, based on the sine-mixture
# simulation from demo/testSine.R.
#
# - "cumr" and "cumc" (and "cumcuda" when built with CUDA) must agree
#   with a brute-force cor() reference within 1e-10.
# - For tInc == 1 they must also agree with the legacy "r" method.
#   (For tInc > 1 the legacy "r" method shifts lags by the column index
#   instead of i*tInc, so it is only used as reference for tInc == 1.)
# - NA input must raise an error for the cum* methods.
# ---------------------------------------------------------------------

library(wcc)

set.seed(260224)

seriesLength <- 1000

makeDyad <- function(seriesLength) {
    t1 <- sin(c(1:seriesLength)/runif(1, min=5, max=20)) + rnorm(seriesLength, mean=0, sd=.5)
    t2 <- sin(c(1:seriesLength)/runif(1, min=5, max=20)) + rnorm(seriesLength, mean=0, sd=.5)
    lagoffset <- floor(runif(1, min=5, max=20))
    predictor <- c(1:seriesLength)/runif(1, min=20, max=30)
    s1 <- s2 <- rep(NA_real_, seriesLength)
    s1[1:(lagoffset-1)] <- rnorm(lagoffset-1, mean=0, sd=.5)
    s2[1:(lagoffset-1)] <- rnorm(lagoffset-1, mean=0, sd=.5)
    tsel1 <- 1:(seriesLength-(lagoffset-1))
    tsel2 <- lagoffset:seriesLength
    s1[tsel2] <- t1[tsel2] * sin(predictor[tsel1]) + t2[tsel1] * cos(predictor[tsel1]) + rnorm(seriesLength-(lagoffset-1), mean=0, sd=.1)
    s2[tsel2] <- t2[tsel2] * cos(predictor[tsel1]) + t1[tsel1] * sin(predictor[tsel1]) + rnorm(seriesLength-(lagoffset-1), mean=0, sd=.1)
    list(s1=s1, s2=s2)
}

# Brute-force reference with explicit lag = i*tInc semantics.
wccBrute <- function(x, y, wMax, tMax, wInc, tInc) {
    n <- length(x)
    tStart <- wMax + tMax
    nRow <- floor((n - tStart) / wInc)
    nLagSteps <- floor(tMax / tInc)
    centerCol <- nLagSteps + 1
    out <- matrix(NA_real_, nrow=nRow, ncol=2*nLagSteps+1)
    for (k in 1:nRow) {
        e <- tStart + (k-1)*wInc
        w <- (e-wMax+1):e
        out[k, centerCol] <- cor(x[w], y[w])
        for (i in 1:nLagSteps) {
            L <- i * tInc
            out[k, centerCol+i] <- cor(x[w], y[w-L])
            out[k, centerCol-i] <- cor(y[w], x[w-L])
        }
    }
    out
}

cudaAvailable <- tryCatch({
    wccCalc(rnorm(200), rnorm(200), wMax=20, tMax=20, method="cumcuda")
    TRUE
}, error = function(e) FALSE)
cat("CUDA backend available:", cudaAvailable, "\n")

maxAbsDiff <- function(a, b) max(abs(a - b), na.rm=TRUE)

params <- expand.grid(wMax=c(50, 100), tMax=c(50, 100), wInc=c(1, 3), tInc=c(1, 2))

nDyads <- 3
for (d in 1:nDyads) {
    dyad <- makeDyad(seriesLength)
    for (p in 1:nrow(params)) {
        wMax <- params$wMax[p]; tMax <- params$tMax[p]
        wInc <- params$wInc[p]; tInc <- params$tInc[p]
        ref  <- wccBrute(dyad$s1, dyad$s2, wMax, tMax, wInc, tInc)
        cumr <- wccCalc(dyad$s1, dyad$s2, wMax=wMax, tMax=tMax, wInc=wInc, tInc=tInc, method="cumr")
        cumc <- wccCalc(dyad$s1, dyad$s2, wMax=wMax, tMax=tMax, wInc=wInc, tInc=tInc, method="cumc")
        stopifnot(all(dim(cumr) == dim(ref)), all(dim(cumc) == dim(ref)))
        dr <- maxAbsDiff(cumr, ref)
        dc <- maxAbsDiff(cumc, ref)
        cat(sprintf("dyad %d wMax=%3d tMax=%3d wInc=%d tInc=%d: |cumr-ref|=%.2e |cumc-ref|=%.2e\n",
                    d, wMax, tMax, wInc, tInc, dr, dc))
        stopifnot(dr < 1e-10, dc < 1e-10)
        if (tInc == 1) {
            legacyR <- wccCalc(dyad$s1, dyad$s2, wMax=wMax, tMax=tMax, wInc=wInc, tInc=tInc, method="r")
            stopifnot(maxAbsDiff(cumr, legacyR) < 1e-10)
        }
        if (cudaAvailable) {
            cumcuda <- wccCalc(dyad$s1, dyad$s2, wMax=wMax, tMax=tMax, wInc=wInc, tInc=tInc, method="cumcuda")
            stopifnot(maxAbsDiff(cumcuda, ref) < 1e-10)
            cumcuda32 <- wccCalc(dyad$s1, dyad$s2, wMax=wMax, tMax=tMax, wInc=wInc, tInc=tInc,
                                 method="cumcuda", precision="single")
            stopifnot(maxAbsDiff(cumcuda32, ref) < 1e-4)
        }
    }
}

# NA input must fail for cum* methods.
naSeries <- rnorm(500)
naSeries[100] <- NA
for (m in c("cumr", "cumc", if (cudaAvailable) "cumcuda")) {
    res <- tryCatch({ wccCalc(naSeries, rnorm(500), wMax=20, tMax=20, method=m); "no error" },
                    error = function(e) "error")
    stopifnot(res == "error")
}

# Zero-variance window must give NA, not a spurious correlation.
# Reset seed so that the constant window is exactly `rep(1, 200)` regardless
# of how much randomness the preceding tests consumed (in particular, the
# cudaAvailable probe above).
set.seed(260224)
constSeries <- c(rnorm(100), rep(1, 200), rnorm(200))
zr <- wccCalc(constSeries, rnorm(500), wMax=20, tMax=20, method="cumr")
zc <- wccCalc(constSeries, rnorm(500), wMax=20, tMax=20, method="cumc")
stopifnot(anyNA(zr), anyNA(zc))

# Deprecated windcross argument still works, with a warning.
w <- tryCatch({ wccCalc(rnorm(300), rnorm(300), wMax=20, tMax=20, windcross=FALSE); "no warning" },
              warning = function(w) "warning")
stopifnot(w == "warning")

# ----------------------------------
# Batched backend: each slab must equal the single-dyad result for its pair.

D <- 5
batchLen <- 800
arr1 <- t(sapply(1:D, function(i) makeDyad(batchLen)$s1))
arr2 <- t(sapply(1:D, function(i) makeDyad(batchLen)$s2))
# mixed pair list: real dyads plus surrogate-style cross pairings
testPairs <- rbind(cbind(1:D, 1:D), c(1, 3), c(4, 2), c(5, 1))

for (m in c("cumc", if (cudaAvailable) "cumcuda")) {
    for (p in 1:nrow(params)) {
        wMax <- params$wMax[p]; tMax <- params$tMax[p]
        wInc <- params$wInc[p]; tInc <- params$tInc[p]
        grids <- wccCalcBatch(arr1, arr2, pairs=testPairs, wMax=wMax, tMax=tMax, wInc=wInc, tInc=tInc, method=m)
        stopifnot(dim(grids)[3] == nrow(testPairs))
        for (q in 1:nrow(testPairs)) {
            single <- wccCalc(arr1[testPairs[q,1],], arr2[testPairs[q,2],],
                              wMax=wMax, tMax=tMax, wInc=wInc, tInc=tInc, method="cumr")
            stopifnot(maxAbsDiff(grids[,,q], single) < 1e-10)
        }
    }
    cat("Batch backend", m, "matches single-dyad results.\n")
}

# FP32 batch path: within 1e-4 of the FP64 reference.
if (cudaAvailable) {
    grids64 <- wccCalcBatch(arr1, arr2, pairs=testPairs, wMax=50, tMax=50, method="cumcuda")
    grids32 <- wccCalcBatch(arr1, arr2, pairs=testPairs, wMax=50, tMax=50,
                            method="cumcuda", precision="single")
    d32 <- maxAbsDiff(grids32, grids64)
    cat(sprintf("FP32 vs FP64 batch: max |diff| = %.2e\n", d32))
    stopifnot(d32 < 1e-4)
}

# precision="single" must be rejected for non-CUDA methods.
res <- tryCatch({ wccCalcBatch(arr1, arr2, wMax=20, tMax=20, method="cumc", precision="single"); "no error" },
                error = function(e) "error")
stopifnot(res == "error")
res <- tryCatch({ wccCalc(rnorm(300), rnorm(300), wMax=20, tMax=20, method="cumc", precision="single"); "no error" },
                error = function(e) "error")
stopifnot(res == "error")

# Default pairs argument: one slab per row, real dyads.
gridsDefault <- wccCalcBatch(arr1, arr2, wMax=50, tMax=50)
stopifnot(dim(gridsDefault)[3] == D,
          maxAbsDiff(gridsDefault[,,2], wccCalc(arr1[2,], arr2[2,], wMax=50, tMax=50, method="cumr")) < 1e-10)

# NA input must fail for the batch backends.
arrNA <- arr1
arrNA[2, 17] <- NA
res <- tryCatch({ wccCalcBatch(arrNA, arr2, wMax=20, tMax=20); "no error" }, error = function(e) "error")
stopifnot(res == "error")

# ----------------------------------
# wccFindDyadParam smoke test on the batched path.

suppressMessages(suppressWarnings({
    fdp <- wccFindDyadParam(inArray1=arr1, inArray2=arr2, wMaxvector=c(50), tMaxvector=c(50),
                            wIncvector=c(1), tIncvector=c(1), Lsizevector=c(8), pspanvector=c(.25),
                            type="Max", nSurrogates=4, method="cumc")
}))
stopifnot(nrow(fdp) == 1,
          all(c("maxMean", "totalMeanKS", "zeroLagVarQdiff") %in% names(fdp)),
          is.finite(fdp$totalMean), is.finite(fdp$maxMeanQdiff))
cat("wccFindDyadParam batched path smoke test passed.\n")

cat("All wccCalc method tests passed.\n")
