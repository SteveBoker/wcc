#
#   Copyright 2001-2026 by the individuals mentioned in the source code history
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

# ---------------------------------------------------------------------
# Program: wccCalcPeakBatch.R
#
# Fused WCC + peak-pick pipeline.  Equivalent to
#
#     g <- wccCalcBatch(seriesArray1, seriesArray2, pairs, ...)
#     g[is.na(g)] <- 0
#     wccPeakPickBatch(g, Lsize, pspan, type, ...)
#
# but with method="cumcuda" the WCC grid never leaves the GPU: the grid
# is computed on the device (zero-variance cells written as 0), smoothed
# with one strided-batched GEMM, searched by the peak kernel, and only
# the nRow x P index/value matrices are downloaded.  This avoids the
# grid's device -> host -> device round trip of the two-step pipeline.
#
# Returns list(index, value): two nRow x P matrices.
# ---------------------------------------------------------------------

wccCalcPeakBatch <- function(seriesArray1, seriesArray2, pairs=NULL,
                             wMax=50, tMax=50, wInc=1, tInc=1,
                             Lsize=8, pspan=.25, type="Max",
                             method=c("cumcuda", "cumc"),
                             precision=c("double", "single")) {
    method <- match.arg(method)
    precision <- match.arg(precision)
    if (precision == "single" && method != "cumcuda") {
        stop(paste0("Warning: precision=\"single\" is only supported with method=\"cumcuda\"."))
    }
    if (!is.numeric(seriesArray1) | !is.numeric(seriesArray2) | !is.matrix(seriesArray1) | !is.matrix(seriesArray2)) {
        stop(paste0("Warning: seriesArray1 and seriesArray2 must be numeric matrices."))
    }
    if (nrow(seriesArray1) != nrow(seriesArray2) && is.null(pairs)) {
        stop(paste0("Warning: seriesArray1 and seriesArray2 must have the same number of rows when pairs is not given."))
    }
    if (ncol(seriesArray1) != ncol(seriesArray2)) {
        stop(paste0("Warning: seriesArray1 and seriesArray2 must have the same number of columns."))
    }
    if (!is.numeric(wMax) | wMax < 5) {
        stop(paste0("Warning: wMax must be a numeric greater than or equal to 5."))
    }
    if (!is.numeric(wInc) | wInc < 1) {
        stop(paste0("Warning: wInc must be a numeric greater than or equal to 1."))
    }
    if (anyNA(seriesArray1) || anyNA(seriesArray2)) {
        stop(paste0("Warning: wccCalcPeakBatch does not support missing data."))
    }
    if (type != "Max" && type != "Min" && type != "max" && type != "min") {
        stop("valid types are: max|Max or Min|min \n")
    }
    if (pspan < 0 || pspan > 1) {
        stop("pspan should be >0 and <1\n")
    }
    if (is.null(pairs)) {
        pairs <- cbind(1:nrow(seriesArray1), 1:nrow(seriesArray1))
    }
    if (!is.matrix(pairs) || ncol(pairs) != 2) {
        stop(paste0("Warning: pairs must be a matrix with two columns."))
    }
    storage.mode(pairs) <- "integer"
    if (anyNA(pairs) || any(pairs[,1] < 1) || any(pairs[,1] > nrow(seriesArray1)) ||
        any(pairs[,2] < 1) || any(pairs[,2] > nrow(seriesArray2))) {
        stop(paste0("Warning: pairs contains an index outside the series arrays."))
    }
    findMin <- (type == "Min" || type == "min")

    if (method == "cumc") {
        g <- wccCalcBatch(seriesArray1, seriesArray2, pairs=pairs,
                          wMax=wMax, tMax=tMax, wInc=wInc, tInc=tInc,
                          method="cumc")
        g[is.na(g)] <- 0
        return(wccPeakPickBatch(g, Lsize=Lsize, pspan=pspan, type=type,
                                method="cumc"))
    }

    storage.mode(seriesArray1) <- "double"
    storage.mode(seriesArray2) <- "double"
    nCol <- 2 * (tMax %/% tInc) + 1
    tLsize <- floor((1/2) * nCol)
    if (Lsize < 1 || Lsize > tLsize) {
        stop(paste("Lsize should be >0 and <= ", tLsize, sep=""))
    }
    M <- wccSmoothMatrix(nCol, pspan)
    res <- .Call("wccpipeline_cuda_batch", seriesArray1, seriesArray2, pairs,
                 as.numeric(wMax), as.numeric(tMax), as.numeric(wInc), as.numeric(tInc),
                 M, as.numeric(Lsize), as.numeric(!findMin),
                 as.integer(precision == "single"),
                 PACKAGE = "wcc")
    P <- nrow(pairs)
    list(index=matrix(res[[1]], ncol=P),
         value=matrix(res[[2]], ncol=P))
}
