# ---------------------------------------------------------------------
# Program: testWccPeakPick.R
#
# Verifies that the matrix-smoother rewrite of wccPeakPick reproduces
# the original per-row loess + spline implementation exactly (up to
# floating point), on WCC grids from the testSine simulation.
# ---------------------------------------------------------------------

library(wcc)

set.seed(260224)

# The original implementation, embedded verbatim as the reference
# (plotting code removed; it was dead, graphs was hard-coded to 0).
wccPeakPickLegacy <- function(tAllCor=NA, Lsize=8, pspan=.25, type="Max") {
    colLen <- length(tAllCor[1,])
    rowLen <- length(tAllCor[,1])
    mx <- rep(NA, (2*colLen-1))
    tIndex <- rep(NA, rowLen)
    tValue <- rep(NA, rowLen)
    findMin <- (type=="Min" || type=="min")
    for (rowNo in c(1:rowLen)) {
        if (any(is.na(tAllCor[rowNo, ]))) next
        tCor <- tAllCor[rowNo, ]
        t1 <- predict(loess(tCor~c(1:colLen), degree=2, span=pspan))
        t2 <- spline(c(1:colLen), t1, n=(2*colLen-1))$y
        windowWidth <- 0
        lookAhead <- 0
        for (j in 1:(colLen-1)) {
            windowWidth <- windowWidth + 1
            tSelect <- (colLen - windowWidth):(colLen + windowWidth)
            mx[j] <- if (findMin) min(t2[tSelect], na.rm=TRUE) else max(t2[tSelect], na.rm=TRUE)
            if (j == 1) mmx <- mx[j]
            else {
                better <- if (findMin) (mx[j] < mmx) else (mx[j] > mmx)
                if (better) {
                    lookAhead <- 0
                    mmx <- mx[j]
                }
                else {
                    lookAhead <- lookAhead + 1
                    if (lookAhead >= Lsize) break
                }
            }
        }
        Index <- match(mmx, t2[tSelect]) + tSelect[1] - 1
        position <- Index - colLen
        if (position > (colLen - Lsize - 1) || position < (-(colLen - Lsize - 1))) {
            tIndex[rowNo] <- NA
            tValue[rowNo] <- NA
        }
        else {
            tIndex[rowNo] <- position
            tValue[rowNo] <- mmx
        }
    }
    if (findMin) return(list(minIndex=tIndex, minValue=tValue))
    return(list(maxIndex=tIndex, maxValue=tValue))
}

makeSeries <- function(n) {
    sin(c(1:n)/runif(1, min=5, max=20)) + rnorm(n, mean=0, sd=.5)
}

sameOrBothNA <- function(a, b, tol=1e-10) {
    all(is.na(a) == is.na(b)) && all(abs(a - b) <= tol, na.rm=TRUE)
}

n <- 1000
for (d in 1:2) {
    s1 <- makeSeries(n)
    s2 <- makeSeries(n)
    for (wMax in c(50, 100)) {
        for (tMax in c(50, 100)) {
            grid <- wccCalc(s1, s2, wMax=wMax, tMax=tMax, method="cumc")
            grid[is.na(grid)] <- 0
            for (type in c("Max", "Min")) {
                for (pspan in c(.25, .4)) {
                    new <- wccPeakPick(grid, Lsize=8, pspan=pspan, type=type)
                    old <- wccPeakPickLegacy(grid, Lsize=8, pspan=pspan, type=type)
                    stopifnot(sameOrBothNA(new[[1]], old[[1]], tol=0),
                              sameOrBothNA(new[[2]], old[[2]]))
                }
            }
            cat(sprintf("dyad %d wMax=%3d tMax=%3d: new == legacy\n", d, wMax, tMax))
        }
    }
}

# Rows containing NA must yield NA index and value.
gridNA <- wccCalc(makeSeries(n), makeSeries(n), wMax=50, tMax=50, method="cumc")
gridNA[is.na(gridNA)] <- 0
gridNA[5, 17] <- NA
res <- wccPeakPick(gridNA, Lsize=8, pspan=.25, type="Max")
stopifnot(is.na(res$maxIndex[5]), is.na(res$maxValue[5]),
          !is.na(res$maxIndex[6]))

# Cache: second build of the same smoother matrix must come from the cache.
key <- paste0(ncol(gridNA), "_", 0.25)
stopifnot(!is.null(wcc:::.wccPeakPickCache[[key]]))
t1 <- system.time(wcc:::wccSmoothMatrix(ncol(gridNA), 0.25))["elapsed"]
stopifnot(t1 < 0.1)

# wccPeakPickBatch (CPU path) must match per-slab wccPeakPick.
arr1 <- t(sapply(1:3, function(i) makeSeries(n)))
arr2 <- t(sapply(1:3, function(i) makeSeries(n)))
grids <- wccCalcBatch(arr1, arr2, wMax=50, tMax=50, method="cumc")
grids[is.na(grids)] <- 0
bp <- wccPeakPickBatch(grids, Lsize=8, pspan=.25, type="Max", method="cumc")
for (p in 1:3) {
    single <- wccPeakPick(grids[,,p], Lsize=8, pspan=.25, type="Max")
    stopifnot(sameOrBothNA(bp$index[,p], single$maxIndex, tol=0),
              sameOrBothNA(bp$value[,p], single$maxValue, tol=0))
}

# CUDA batch paths (FP64 and FP32) against the CPU batch result.
cudaAvailable <- tryCatch({
    wccCalc(rnorm(200), rnorm(200), wMax=20, tMax=20, method="cumcuda")
    TRUE
}, error = function(e) FALSE)

if (cudaAvailable) {
    bg <- wccPeakPickBatch(grids, Lsize=8, pspan=.25, type="Max", method="cumcuda")
    stopifnot(sameOrBothNA(bg$index, bp$index, tol=0),
              sameOrBothNA(bg$value, bp$value))

    b32 <- wccPeakPickBatch(grids, Lsize=8, pspan=.25, type="Max",
                            method="cumcuda", precision="single")
    stopifnot(all(is.na(b32$index) == is.na(bp$index)))
    # Peak values within FP32 tolerance everywhere.
    stopifnot(all(abs(b32$value - bp$value) < 1e-4, na.rm=TRUE))
    # Indices match except possibly on near-ties, where a shift of one
    # half-step is acceptable.
    idxDiff <- abs(b32$index - bp$index)
    stopifnot(all(idxDiff <= 1, na.rm=TRUE))
    nShift <- sum(idxDiff > 0, na.rm=TRUE)
    cat(sprintf("FP32 peak pick: %d/%d indices shifted by one half-step\n",
                nShift, sum(!is.na(idxDiff))))
}

# precision="single" must be rejected for method="cumc".
res <- tryCatch({ wccPeakPickBatch(grids, method="cumc", precision="single"); "no error" },
                error = function(e) "error")
stopifnot(res == "error")

# ----------------------------------
# Fused pipeline: wccCalcPeakBatch must match the two-step pipeline.

fpc <- wccCalcPeakBatch(arr1, arr2, wMax=50, tMax=50, Lsize=8, pspan=.25,
                        type="Max", method="cumc")
stopifnot(sameOrBothNA(fpc$index, bp$index, tol=0),
          sameOrBothNA(fpc$value, bp$value, tol=0))

if (cudaAvailable) {
    fg <- wccCalcPeakBatch(arr1, arr2, wMax=50, tMax=50, Lsize=8, pspan=.25,
                           type="Max", method="cumcuda")
    stopifnot(all(dim(fg$index) == dim(bp$index)),
              sameOrBothNA(fg$index, bp$index, tol=0),
              sameOrBothNA(fg$value, bp$value))

    # Repeated calls must keep working (persistent buffers/handle reuse).
    fg2 <- wccCalcPeakBatch(arr1, arr2, wMax=50, tMax=50, Lsize=8, pspan=.25,
                            type="Max", method="cumcuda")
    stopifnot(identical(fg, fg2))

    # A smaller follow-up problem must not read stale buffer tails.
    fgSmall <- wccCalcPeakBatch(arr1[1:2, 1:600], arr2[1:2, 1:600],
                                wMax=30, tMax=30, Lsize=8, pspan=.25,
                                type="Max", method="cumcuda")
    fgSmallRef <- wccCalcPeakBatch(arr1[1:2, 1:600], arr2[1:2, 1:600],
                                   wMax=30, tMax=30, Lsize=8, pspan=.25,
                                   type="Max", method="cumc")
    stopifnot(sameOrBothNA(fgSmall$index, fgSmallRef$index, tol=0),
              sameOrBothNA(fgSmall$value, fgSmallRef$value))

    # FP32 fused path within tolerance.
    f32 <- wccCalcPeakBatch(arr1, arr2, wMax=50, tMax=50, Lsize=8, pspan=.25,
                            type="Max", method="cumcuda", precision="single")
    stopifnot(all(abs(f32$value - bp$value) < 1e-4, na.rm=TRUE),
              all(abs(f32$index - bp$index) <= 1, na.rm=TRUE))

    # "Min" must agree with the two-step CPU pipeline too.
    bpMin <- wccPeakPickBatch(grids, Lsize=8, pspan=.25, type="Min", method="cumc")
    fgMin <- wccCalcPeakBatch(arr1, arr2, wMax=50, tMax=50, Lsize=8, pspan=.25,
                              type="Min", method="cumcuda")
    stopifnot(sameOrBothNA(fgMin$index, bpMin$index, tol=0),
              sameOrBothNA(fgMin$value, bpMin$value))
    cat("Fused wccCalcPeakBatch (cumcuda) matches the two-step pipeline.\n")
}

# precision="single" must be rejected for method="cumc" in the fused path.
res <- tryCatch({ wccCalcPeakBatch(arr1, arr2, method="cumc", precision="single"); "no error" },
                error = function(e) "error")
stopifnot(res == "error")

cat("All wccPeakPick tests passed.\n")
