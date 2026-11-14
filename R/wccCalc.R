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
# Program: wccCalc.R
#  Author: Steven Boker
#    Date: Wed Jan 28 08:52:17 EST 2026
#
# This function runs a windowed cross correlation and returns a matrix whose columns
# are the lags and whose rows are the elapsed time of each window.
#
# Backends are selected with the method argument:
#   "cumr"    -- cumulative-sum algorithm, pure R (default)
#   "cumc"    -- cumulative-sum algorithm, C
#   "cumcuda" -- cumulative-sum algorithm, CUDA (requires a CUDA-enabled build)
#   "c"       -- original C implementation (windcross)
#   "r"       -- original pure-R loop
#
# The cum* methods compute each cell in O(1) from prefix sums, so their
# cost is O(n * nLags) independent of window size.  They do not support
# missing data: any NA in the input is an error.
#
# Note on lag stepping when tInc > 1: the cum* methods (like "c") use
# lags 0, tInc, 2*tInc, ..., tMax.  The legacy "r" method shifts by the
# column index itself (lags 0, 1, 2, ...) regardless of tInc.
#
# ---------------------------------------------------------------------
# Revision History
#  Steve Boker  -- Wed Jan 28 08:52:19 EST 2026
#      Created wccCalc.R
#
# ---------------------------------------------------------------------

# ----------------------------------
# Calculate WCC.

wccCalc <- function(inSeries1, inSeries2, wMax=50, tMax=50, wInc=1, tInc=1,
                    method=c("c", "cumr", "cumc", "cumcuda", "r"),
                    precision=c("double", "single"), ...) {
    # Deprecation: allow old windcross argument
    dots <- list(...)
    if ("windcross" %in% names(dots)) {
        warning("Argument 'windcross' is deprecated; use method=\"c\" or method=\"r\" instead.")
        method <- if (isTRUE(dots$windcross)) "c" else "r"
    } else {
        method <- match.arg(method)
    }
    precision <- match.arg(precision)
    if (precision == "single" && method != "cumcuda") {
        stop(paste0("Warning: precision=\"single\" is only supported with method=\"cumcuda\"."))
    }
    if (!is.numeric(inSeries1) | !is.numeric(inSeries2) | !is.vector(inSeries1) | !is.vector(inSeries2) | length(inSeries1) != length(inSeries2)) {
        stop(paste0("Warning: inSeries1 and inSeries2 must be numeric vectors of equal length."))
    }
    if (!is.numeric(wMax) | wMax < 5  ) {
        stop(paste0("Warning: wMax must be a numeric greater than or equal to 5."))
    }
    if ( !is.numeric(wInc) |  wInc < 1 ) {
        stop(paste0("Warning: wInc must be a numeric greater than or equal to 1."))
    }

    maxRowLen <- length(inSeries1)
    tStart <- wMax + tMax
	nRow <- floor((maxRowLen - tStart) / wInc)
	nCol <- 2*(floor(tMax / tInc)) + 1
    centerCol <- (floor(tMax / tInc)) + 1
    if (nRow < 0) {
        stop(paste0("Warning: bad choice for wMax and/or wInc parameters. The result matrix has ", nRow, " rows."))
    }
    if (nCol < 0) {
        stop(paste0("Warning: bad choice for tMax and/or tInc  parameters. The result matrix has ", nCol, " columns."))
    }
    if (method %in% c("cumr", "cumc", "cumcuda")) {
        if (anyNA(inSeries1) || anyNA(inSeries2)) {
            stop(paste0("Warning: method \"", method, "\" does not support missing data. Use method=\"c\" or method=\"r\"."))
        }
    }
    if (method == "cumr") {
        return(wccCalcCumR(as.numeric(inSeries1), as.numeric(inSeries2), wMax, tMax, wInc, tInc))
    }
    if (method == "cumc") {
        return(.Call("windcrosscum", as.numeric(inSeries1), as.numeric(inSeries2), as.numeric(wMax), as.numeric(tMax), as.numeric(wInc), as.numeric(tInc), PACKAGE = "wcc"))
    }
    if (method == "cumcuda") {
        return(.Call("windcrosscum_cuda", as.numeric(inSeries1), as.numeric(inSeries2), as.numeric(wMax), as.numeric(tMax), as.numeric(wInc), as.numeric(tInc), as.integer(precision == "single"), PACKAGE = "wcc"))
    }
    if (method == "c") {
        tData <- .Call("windcross", as.numeric(inSeries1), as.numeric(inSeries2), as.numeric(wMax), as.numeric(tMax), as.numeric(wInc), as.numeric(tInc), PACKAGE = "wcc")
        tData[tData < -99] <- NA
        return(tData)
    }
    else {
        tData <- matrix(NA, nrow=nRow, ncol=nCol)
        elapsedSeq <-  seq(from=tStart, to=maxRowLen, by=wInc)
        for (tRow in 1:nRow) {
            elapsedIndex <- elapsedSeq[tRow]
            lag0Sel <- seq(from=elapsedIndex, to=(1+ elapsedIndex-wMax), by=-1)
            lag0win1 <- inSeries1[lag0Sel]
            lag0win2 <- inSeries2[lag0Sel]
            tData[tRow,centerCol] <- cor(lag0win1,lag0win2, use="pairwise.complete.obs")
            for(i in 1:(floor(tMax / tInc))) {
                tData[tRow,(centerCol+i)] <- cor(lag0win1, inSeries2[lag0Sel-i], use="pairwise.complete.obs")
                tData[tRow,(centerCol-i)] <- cor(lag0win2, inSeries1[lag0Sel-i], use="pairwise.complete.obs")
            }
        }
        return(tData)
    }
}

# ----------------------------------
# Cumulative-sum WCC in pure R.
#
# All window sums are differences of prefix sums, so each cell costs O(1)
# and all rows of a lag column are computed in one vectorized step.
# Series are centered by their global means first to limit floating-point
# cancellation in the prefix-sum differences.

wccCalcCumR <- function(x, y, wMax, tMax, wInc, tInc) {
    n <- length(x)
    tStart <- wMax + tMax
    nRow <- floor((n - tStart) / wInc)
    nLagSteps <- floor(tMax / tInc)
    nCol <- 2 * nLagSteps + 1
    centerCol <- nLagSteps + 1

    x <- x - mean(x)
    y <- y - mean(y)

    # Prefix sums with a leading zero so that
    # sum(v[a:b]) == cs[b + 1] - cs[a].
    csX  <- c(0, cumsum(x))
    csY  <- c(0, cumsum(y))
    csX2 <- c(0, cumsum(x * x))
    csY2 <- c(0, cumsum(y * y))

    # Row windows end at hi and start at lo + 1 (window length wMax).
    hi <- tStart + (seq_len(nRow) - 1) * wInc
    lo <- hi - wMax

    # Base-window statistics, hoisted: computed once per row, reused for
    # every lag column.
    bSx <- csX[hi + 1] - csX[lo + 1]
    bQx <- csX2[hi + 1] - csX2[lo + 1]
    bSy <- csY[hi + 1] - csY[lo + 1]
    bQy <- csY2[hi + 1] - csY2[lo + 1]
    bVarx <- wMax * bQx - bSx * bSx
    bVary <- wMax * bQy - bSy * bSy

    tData <- matrix(NA_real_, nrow = nRow, ncol = nCol)

    # Center column: lag 0.
    csXY <- c(0, cumsum(x * y))
    sxy <- csXY[hi + 1] - csXY[lo + 1]
    num <- wMax * sxy - bSx * bSy
    den <- bVarx * bVary
    tData[, centerCol] <- ifelse(den > 0, num / sqrt(den), NA_real_)

    if (nLagSteps > 0) {
        for (i in seq_len(nLagSteps)) {
            L <- i * tInc
            # Lagged-window prefix indices: window [lo+1-L, hi-L].
            hiL <- hi - L
            loL <- lo - L

            # Column centerCol + i: cor(x[window], y[window - L]).
            p <- x * c(rep(0, L), y[seq_len(n - L)])
            csP <- c(0, cumsum(p))
            sxy <- csP[hi + 1] - csP[lo + 1]
            sy  <- csY[hiL + 1] - csY[loL + 1]
            qy  <- csY2[hiL + 1] - csY2[loL + 1]
            vary <- wMax * qy - sy * sy
            num <- wMax * sxy - bSx * sy
            den <- bVarx * vary
            tData[, centerCol + i] <- ifelse(den > 0, num / sqrt(den), NA_real_)

            # Column centerCol - i: cor(y[window], x[window - L]).
            p <- y * c(rep(0, L), x[seq_len(n - L)])
            csP <- c(0, cumsum(p))
            sxy <- csP[hi + 1] - csP[lo + 1]
            sx  <- csX[hiL + 1] - csX[loL + 1]
            qx  <- csX2[hiL + 1] - csX2[loL + 1]
            varx <- wMax * qx - sx * sx
            num <- wMax * sxy - bSy * sx
            den <- varx * bVary
            tData[, centerCol - i] <- ifelse(den > 0, num / sqrt(den), NA_real_)
        }
    }

    tData
}
