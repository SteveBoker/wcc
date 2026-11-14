# ---------------------------------------------------------------------
# Program: benchmarkWCC.R
#
# Benchmarks the wccCalc backends ("cumr", "cumc", "cumcuda", "c", "r")
# on sine-mixture series as in testSine.R, and reports timings and
# speedups relative to the original C implementation ("c").
#
# The slow legacy methods ("c", "r") are skipped for the largest series.
# ---------------------------------------------------------------------

library(wcc)

set.seed(260224)

makeSeries <- function(n) {
    sin(c(1:n)/runif(1, min=5, max=20)) + rnorm(n, mean=0, sd=.5)
}

cudaAvailable <- tryCatch({
    wccCalc(rnorm(200), rnorm(200), wMax=20, tMax=20, method="cumcuda")
    TRUE
}, error = function(e) FALSE)

methods <- c("cumr", "cumc", if (cudaAvailable) "cumcuda", "c", "r")

lengths <- c(1e3, 1e4, 1e5)
params <- expand.grid(wMax=c(50, 100), tMax=c(50, 100))

results <- data.frame()
for (n in lengths) {
    s1 <- makeSeries(n)
    s2 <- makeSeries(n)
    for (p in 1:nrow(params)) {
        wMax <- params$wMax[p]; tMax <- params$tMax[p]
        row <- data.frame(n=n, wMax=wMax, tMax=tMax)
        for (m in methods) {
            # legacy methods are O(nRow * nLag * wMax); cap their runtime
            if (m %in% c("c", "r") && n > 1e4) {
                row[[m]] <- NA_real_
                next
            }
            tm <- system.time(wccCalc(s1, s2, wMax=wMax, tMax=tMax, method=m))["elapsed"]
            row[[m]] <- as.numeric(tm)
        }
        results <- rbind(results, row)
        cat(sprintf("n=%6d wMax=%3d tMax=%3d : %s\n", n, wMax, tMax,
                    paste(sprintf("%s=%.3fs", methods, unlist(row[methods])), collapse="  ")))
    }
}

cat("\nSpeedup relative to method=\"c\" (where measured):\n")
speedup <- results
for (m in methods) speedup[[m]] <- results$c / results[[m]]
print(speedup, digits=3)

# ----------------------------------
# Batched backend vs a loop of single-dyad calls.

batchMethods <- c("cumc", if (cudaAvailable) "cumcuda")
cat("\nBatched wccCalcBatch vs single-dyad loop (wMax=tMax=50):\n")
for (D in c(100, 500)) {
    for (n in c(1000, 5000)) {
        arr1 <- t(sapply(1:D, function(i) makeSeries(n)))
        arr2 <- t(sapply(1:D, function(i) makeSeries(n)))
        tLoop <- system.time(for (i in 1:D) wccCalc(arr1[i,], arr2[i,], wMax=50, tMax=50, method="cumc"))["elapsed"]
        for (m in batchMethods) {
            tBatch <- system.time(wccCalcBatch(arr1, arr2, wMax=50, tMax=50, method=m))["elapsed"]
            cat(sprintf("D=%4d n=%5d %-8s: batch=%.2fs  loop(cumc)=%.2fs  speedup=%.1fx\n",
                        D, n, m, tBatch, tLoop, tLoop / tBatch))
        }
    }
}

# ----------------------------------
# Peak picking: matrix smoother vs the original per-row loess.

ppGrid <- wccCalc(makeSeries(10000), makeSeries(10000), wMax=100, tMax=100, method="cumc")
ppGrid[is.na(ppGrid)] <- 0
invisible(wccPeakPick(ppGrid[1:2,], Lsize=8, pspan=.25))  # warm the smoother cache
tNewPP <- system.time(wccPeakPick(ppGrid, Lsize=8, pspan=.25))["elapsed"]
tOldPP <- system.time({
    for (r in 1:nrow(ppGrid)) {
        t1 <- predict(loess(ppGrid[r,] ~ c(1:ncol(ppGrid)), degree=2, span=.25))
        t2 <- spline(c(1:ncol(ppGrid)), t1, n=(2*ncol(ppGrid)-1))$y
    }
})["elapsed"]
cat(sprintf("\nwccPeakPick on a %d x %d grid: matrix smoother=%.2fs  per-row loess (smoothing only)=%.2fs  speedup=%.0fx\n",
            nrow(ppGrid), ncol(ppGrid), tNewPP, tOldPP, tOldPP / tNewPP))

# ----------------------------------
# wccFindDyadParam: batched path vs legacy per-dyad path.

D <- 20
n <- 1000
arr1 <- t(sapply(1:D, function(i) makeSeries(n)))
arr2 <- t(sapply(1:D, function(i) makeSeries(n)))
suppressMessages({
    tNew <- system.time(wccFindDyadParam(inArray1=arr1, inArray2=arr2, wMaxvector=c(50), tMaxvector=c(50),
                        nSurrogates=10, method="cumc"))["elapsed"]
    tOld <- system.time(wccFindDyadParam(inArray1=arr1, inArray2=arr2, wMaxvector=c(50), tMaxvector=c(50),
                        nSurrogates=10, method="r"))["elapsed"]
})
cat(sprintf("\nwccFindDyadParam (D=%d, n=%d, 10 surrogates): batched cumc=%.1fs  legacy r=%.1fs  speedup=%.1fx\n",
            D, n, tNew, tOld, tOld / tNew))
