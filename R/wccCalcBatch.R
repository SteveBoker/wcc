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
# Program: wccCalcBatch.R
#
# Batched windowed cross correlation: computes the WCC grid for a list
# of series pairs in one call.  Per-series prefix sums are computed once
# and shared across all pairs, so real dyads and surrogate pairings of
# the same series pool can be evaluated together cheaply.
#
# seriesArray1, seriesArray2: D x n numeric matrices (rows = series).
# pairs: P x 2 integer matrix; row p selects series pairs[p,1] from
#        seriesArray1 and pairs[p,2] from seriesArray2.  Defaults to the
#        D real dyads cbind(1:D, 1:D).
#
# Returns an nRow x nCol x P array; slab [,,p] equals
# wccCalc(seriesArray1[pairs[p,1],], seriesArray2[pairs[p,2],], ...).
#
# Missing data is not supported: any NA in the inputs is an error.
# ---------------------------------------------------------------------

wccCalcBatch <- function(seriesArray1, seriesArray2, pairs=NULL,
                         wMax=50, tMax=50, wInc=1, tInc=1,
                         method=c("cumc", "cumcuda"),
                         precision=c("double", "single")) {
    # Note: wccCalcBatch only supports batched backends (cumc/cumcuda),
    # not the single-dyad paths. Use method="cumc" for the default batched C backend.
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
    if (!is.numeric(wMax) | wMax < 5  ) {
        stop(paste0("Warning: wMax must be a numeric greater than or equal to 5."))
    }
    if ( !is.numeric(wInc) |  wInc < 1 ) {
        stop(paste0("Warning: wInc must be a numeric greater than or equal to 1."))
    }
    if (anyNA(seriesArray1) || anyNA(seriesArray2)) {
        stop(paste0("Warning: wccCalcBatch does not support missing data."))
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
    storage.mode(seriesArray1) <- "double"
    storage.mode(seriesArray2) <- "double"

    if (method == "cumc") {
        .Call("windcrosscum_batch", seriesArray1, seriesArray2, pairs,
              as.numeric(wMax), as.numeric(tMax), as.numeric(wInc), as.numeric(tInc),
              PACKAGE = "wcc")
    } else {
        .Call("windcrosscum_cuda_batch", seriesArray1, seriesArray2, pairs,
              as.numeric(wMax), as.numeric(tMax), as.numeric(wInc), as.numeric(tInc),
              as.integer(precision == "single"),
              PACKAGE = "wcc")
    }
}
