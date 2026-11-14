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
# Program: peakPick.R
#  Author: Steven Boker, Minquan Xu
#    Date: Wed Jan 28 08:55:03 EST 2026
#
# Finds, for each row of a WCC matrix, the local peak (or valley)
# closest to zero lag.
#
# Each row is smoothed with loess(degree=2, span=pspan) evaluated at the
# observed lags and then interpolated to 2*colLen-1 points with a cubic
# spline.  For fixed predictor positions both operations are linear in
# the row values, so the combined smoother is a single
# (2*colLen-1) x colLen matrix that depends only on (colLen, pspan).
# wccSmoothMatrix builds that matrix once (by pushing the columns of the
# identity matrix through the original loess+spline pipeline, which
# reproduces it exactly) and caches it; a whole WCC matrix is then
# smoothed with one matrix product instead of one loess fit per row.
#
# ---------------------------------------------------------------------
# Revision History
#  Steve Boker  -- 2001
#      Created peakPick.R
#
# ---------------------------------------------------------------------

# input data structure:
#    row: time
#    column: cross correlations
# output: a list of local peak indice and values

# Parameters:
# ---------------------------------------------------------------------
# tAllCor: numeric matrix of windowed cross correlations (rows = windows,
#          columns = lags)
# Lsize: local search region, the value should be larger than 0 and
#        less than 1/2 length of one row
# pspan: see loess() for span
# type: local maximum or local minimum, valid values are: "Min" and "Max"
#---------------------------------------------------------------------

# Cache of smoother matrices, keyed by "colLen_pspan".
.wccPeakPickCache <- new.env(parent=emptyenv())

# ----------------------------------
# The combined loess + spline smoother as an explicit linear operator.
# Returns the (2*colLen-1) x colLen matrix M such that
# M %*% y == spline(1:colLen, predict(loess(y ~ 1:colLen, degree=2,
#                   span=pspan)), n=2*colLen-1)$y  for any y.

wccSmoothMatrix <- function(colLen, pspan) {
    key <- paste0(colLen, "_", pspan)
    M <- .wccPeakPickCache[[key]]
    if (!is.null(M)) {
        return(M)
    }
    M <- matrix(0, nrow=2*colLen-1, ncol=colLen)
    x <- c(1:colLen)
    for (j in 1:colLen) {
        e <- numeric(colLen)
        e[j] <- 1
        t1 <- predict(loess(e ~ x, degree=2, span=pspan))
        M[, j] <- spline(x, t1, n=(2*colLen-1))$y
    }
    .wccPeakPickCache[[key]] <- M
    M
}

# ----------------------------------
# Expanding-window peak search on one smoothed row t2 (length 2*colLen-1,
# center at index colLen).  Replicates the original look-ahead loop:
# mx[j] is the running max over the window of half-width j; the search
# stops after Lsize consecutive non-improvements.  Returns c(index, value)
# with NAs when the peak fails the boundary criterion.

wccPeakSearch <- function(t2, colLen, Lsize) {
    # mx[j] = max(t2[(colLen-j):(colLen+j)]) for all j at once
    lv <- cummax(t2[colLen:1])
    rv <- cummax(t2[colLen:(2*colLen-1)])
    mxAll <- pmax(lv, rv)
    mx <- mxAll[2:colLen]            # mx[j], j = 1 .. colLen-1

    mmx <- mx[1]
    lookAhead <- 0
    windowWidth <- colLen - 1
    if (colLen > 2) {
        for (j in 2:(colLen-1)) {
            if (mx[j] > mmx) {
                lookAhead <- 0
                mmx <- mx[j]
            }
            else {
                lookAhead <- lookAhead + 1
                if (lookAhead >= Lsize) {
                    windowWidth <- j
                    break
                }
            }
        }
    }
    tSelect <- (colLen - windowWidth):(colLen + windowWidth)
    Index <- match(mmx, t2[tSelect]) + tSelect[1] - 1
    position <- Index - colLen
    if (position > (colLen - Lsize - 1) || position < (-(colLen - Lsize - 1))) {
        c(NA_real_, NA_real_)
    }
    else {
        c(position, mmx)
    }
}

wccPeakPick <- function(tAllCor=NA, Lsize=8, pspan=.25, type="Max") {
    #----------------check for validness of parameters -------------------
    colLen <- length(tAllCor[1,]) # col length --- number of columns
    rowLen <- length(tAllCor[,1]) # row length --- number of rows
    tLsize <- floor((1/2)*colLen) # maximum local search region

    if(Lsize<1 || Lsize>tLsize) { # Lsize too small or too large
          errorStr<- paste("Lsize should be >0 and <= ", tLsize, sep="")
          stop(errorStr)
    }
    if(pspan<0 || pspan>1) { # invalid pspan value
           stop("pspan should be >0 and <1\n")
    }
    if(type!="Max"&&type!="Min"&&type!="max"&&type!="min"){ # only two types
           stop("valid types are: max|Max or Min|min \n")
    }
    findMin <- (type=="Min" || type=="min")

    tIndex <- rep(NA, rowLen) #vector of peak index---one peak index for each row
    tValue <- rep(NA, rowLen) #vector of peak value---one peak value for each row

    # Rows containing any NA are skipped, as in the original row loop.
    okRows <- which(rowSums(is.na(tAllCor)) == 0)
    if (length(okRows) > 0) {
        M <- wccSmoothMatrix(colLen, pspan)
        # smooth all NA-free rows with one matrix product
        T2 <- tAllCor[okRows, , drop=FALSE] %*% t(M)
        for (r in seq_along(okRows)) {
            t2 <- T2[r, ]
            if (findMin) {
                res <- wccPeakSearch(-t2, colLen, Lsize)
                tIndex[okRows[r]] <- res[1]
                tValue[okRows[r]] <- -res[2]
            }
            else {
                res <- wccPeakSearch(t2, colLen, Lsize)
                tIndex[okRows[r]] <- res[1]
                tValue[okRows[r]] <- res[2]
            }
        }
    }

    if (findMin) {
        return(list(minIndex=tIndex, minValue=tValue))
    }
    return(list(maxIndex=tIndex, maxValue=tValue))
}

# ----------------------------------
# Batched peak picking over a 3-D array of WCC grids (nRow x nCol x P),
# as returned by wccCalcBatch.  method="cumc" runs the fast CPU path per
# slab; method="cumcuda" smooths all slabs with one batched GEMM on the
# GPU and runs the search kernel there (requires a CUDA-enabled build).
# Returns list(index, value): two nRow x P matrices.

wccPeakPickBatch <- function(grids, Lsize=8, pspan=.25, type="Max", method=c("cumc", "cumcuda"),
                             precision=c("double", "single")) {
    method <- match.arg(method)
    precision <- match.arg(precision)
    if (precision == "single" && method != "cumcuda") {
        stop("precision=\"single\" is only supported with method=\"cumcuda\".")
    }
    if (!is.array(grids) || length(dim(grids)) != 3) {
        stop("grids must be a 3-dimensional array (nRow x nCol x P).")
    }
    nRow <- dim(grids)[1]
    colLen <- dim(grids)[2]
    P <- dim(grids)[3]
    findMin <- (type=="Min" || type=="min")

    if (method == "cumcuda") {
        if (anyNA(grids)) {
            stop("method \"cumcuda\" does not support NA cells; replace them before peak picking.")
        }
        M <- wccSmoothMatrix(colLen, pspan)
        res <- .Call("wccpeakpick_cuda_batch", grids, M,
                     as.numeric(Lsize), as.numeric(!findMin),
                     as.integer(precision == "single"), PACKAGE = "wcc")
        return(list(index=matrix(res[[1]], nrow=nRow, ncol=P),
                    value=matrix(res[[2]], nrow=nRow, ncol=P)))
    }

    index <- matrix(NA_real_, nrow=nRow, ncol=P)
    value <- matrix(NA_real_, nrow=nRow, ncol=P)
    for (p in 1:P) {
        ppOut <- wccPeakPick(grids[,,p], Lsize=Lsize, pspan=pspan, type=type)
        if (findMin) {
            index[,p] <- ppOut$minIndex
            value[,p] <- ppOut$minValue
        }
        else {
            index[,p] <- ppOut$maxIndex
            value[,p] <- ppOut$maxValue
        }
    }
    list(index=index, value=value)
}
