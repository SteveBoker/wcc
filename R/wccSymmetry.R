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
# Program: wccSymmetry.R
#  Author: Steven Boker, Claude
#    Date: Fri Nov  6 2026
#
# Estimates when symmetry forms and breaks in a WCC matrix.
#
# For every cell of the matrix (rows = elapsed time, columns = lags, as
# returned by wccCalc) the direction of zero change -- the local isobar --
# is computed from the gradient (dr/dt, dr/dlag).  Within each row the
# spread of these directions across lags is summarized.  Low spread means
# the whole lag profile is moving the same way: either standing still
# (stationarity) or drifting uniformly in lag.  Both count as symmetry.
#
# Why: if the profile translates rigidly in lag, r(t, lag) = f(lag - v t),
# then dr/dt = -v f' and dr/dlag = f', so the zero-change slope
# dlag/dt = -(dr/dt)/(dr/dlag) = v at every lag.  The directions agree
# across the row and their dispersion is zero.  Stationarity is v = 0.
#
# Directions are axial (a line has no arrow), so angles are doubled before
# averaging.  Each cell is weighted by its squared gradient magnitude since
# direction is undefined where the surface is flat.  With that weighting
# the dispersion is 1 minus the coherence of the row's 2 x 2 structure
# tensor.
#
# Two preprocessing steps matter on real WCC matrices, in this order:
#  1. Standardizing each row across lags so that changes in peak height
#     alone, r = a(t) f(lag), which give lag-dependent zero-change slopes,
#     are not read as broken symmetry.  Done first so low-amplitude windows
#     are not swamped by their neighbors in step 2.
#  2. Smoothing along elapsed time (default wMax / (4 wInc) rows) so
#     window-to-window sampling noise in r is not read as change.  Noisy
#     data may need up to one window (wMax / wInc rows).
#
# wccSymmetryEpochs thresholds the dispersion.  By default the threshold
# is chosen from the data by splitting log dispersion into two groups
# (Otsu's method), because absolute dispersion levels depend strongly on
# the noise in the data.
#
# Gradients are taken in index units (rows x lag columns) so the angles
# do not depend on the arbitrary ratio of time units to lag units;
# wInc, tInc and samplespersecond only rescale the reported times and
# drift rate.
#
# ---------------------------------------------------------------------
# Revision History
#  Steve Boker, Claude  -- Fri Nov  6 2026
#      Created wccSymmetry.R
#      Standardize rows before smoothing; lighter default smoothing;
#      data-driven ("auto") threshold by default.
#
# ---------------------------------------------------------------------

# ----------------------------------
# Per-row dispersion of zero-change directions.

wccZeroChangeDirections <- function(wccMatrix, wMax=NA, tMax=NA, wInc=1,
                                    tInc=1, samplespersecond=1,
                                    method=c("shift", "gradient"),
                                    timeScale=NA, maxShift=NA,
                                    timeSmooth=NA, normalizeRows=TRUE,
                                    weight=c("gradient", "none"),
                                    minGrad=1e-8) {
    method <- match.arg(method)
    weight <- match.arg(weight)
    m <- as.matrix(wccMatrix)
    if (!is.numeric(m)) {
        stop("Warning: wccMatrix must be a numeric matrix as returned by wccCalc.")
    }
    nT <- nrow(m)
    nL <- ncol(m)
    if (nT < 3 || nL < 3) {
        stop("Warning: wccMatrix needs at least 3 rows (windows) and 3 columns (lags).")
    }

    # Row-level descriptions of the raw cross correlations, used to
    # describe epochs (strength, location and sharpness of the peak).
    raw <- m
    raw[!is.finite(raw)] <- -Inf
    am <- max.col(raw, ties.method="first")
    rr <- seq_len(nT)
    peakR <- m[cbind(rr, am)]
    peakSharp <- peakR - 0.5 * (m[cbind(rr, pmax(1, am - 3))] + m[cbind(rr, pmin(nL, am + 3))])
    centerCol <- (nL + 1) / 2
    peakLag <- (am - centerCol) * tInc / samplespersecond

    if (method == "shift") {
        if (is.na(timeScale)) {
            if (is.na(wMax)) {
                stop("Warning: method=\"shift\" needs wMax or timeScale.")
            }
            timeScale <- max(1, round(wMax / (2 * wInc)))
        }
        if (is.na(maxShift)) maxShift <- max(1, ceiling(timeScale / 2))
        if (is.na(timeSmooth)) timeSmooth <- 1
    } else {
        if (is.na(timeSmooth)) {
            timeSmooth <- if (is.na(wMax)) 1 else max(1, round(wMax / (4 * wInc)))
        }
    }

    # 1. Standardize each row across lags.  This is done before smoothing
    #    so that low-amplitude windows (e.g. where one series is noise) are
    #    not swamped by high-amplitude neighbors in the moving average.
    if (normalizeRows) {
        m <- t(apply(m, 1, function(r) {
            s <- stats::sd(r, na.rm=TRUE)
            if (!is.finite(s) || s == 0) return(r * 0)
            (r - mean(r, na.rm=TRUE)) / s
        }))
    }

    # 2. Smooth along elapsed time (centered moving average).
    if (timeSmooth > 1) {
        k <- min(as.integer(timeSmooth), nT)
        m <- apply(m, 2, function(cl) {
            s <- as.numeric(stats::filter(cl, rep(1 / k, k), sides=2))
            s[is.na(s)] <- cl[is.na(s)]
            s
        })
    }

    rowIdx <- seq_len(nT)
    offset <- if (is.na(wMax) || is.na(tMax)) 0 else wMax + tMax
    params <- list(method=method, wMax=wMax, tMax=tMax, wInc=wInc, tInc=tInc,
                   samplespersecond=samplespersecond, timeScale=timeScale)

    if (method == "shift") {
        # Compare the lag profile timeScale windows ahead with the one
        # timeScale windows behind, for every uniform shift in lag up to
        # maxShift columns.  The best shift is the common direction of zero
        # change; 1 - its correlation measures how far the change departs
        # from a uniform translation.
        h <- as.integer(timeScale)
        if (2 * h + 1 > nT) {
            stop("Warning: timeScale is too large for the number of windows.")
        }
        maxShift <- min(as.integer(maxShift), nL - 3)
        idx <- (h + 1):(nT - h)
        best <- rep(-Inf, length(idx))
        bestShift <- rep(0, length(idx))
        for (s in -maxShift:maxShift) {
            a <- if (s >= 0) (1 + s):nL else 1:(nL + s)
            b <- a - s
            A <- m[idx + h, a, drop=FALSE]
            B <- m[idx - h, b, drop=FALSE]
            A <- A - rowMeans(A, na.rm=TRUE)
            B <- B - rowMeans(B, na.rm=TRUE)
            r <- rowSums(A * B, na.rm=TRUE) /
                 sqrt(rowSums(A^2, na.rm=TRUE) * rowSums(B^2, na.rm=TRUE))
            r[!is.finite(r)] <- -1
            better <- r > best
            best[better] <- r[better]
            bestShift[better] <- s
        }
        coh <- rep(NA_real_, nT)
        shiftCols <- rep(NA_real_, nT)
        coh[idx] <- best
        shiftCols[idx] <- bestShift
        dirn <- atan(shiftCols / (2 * h))
        out <- data.frame(
            window = rowIdx,
            elapsedTime = (offset + (rowIdx - 1) * wInc) / samplespersecond,
            dispersion = 1 - coh,
            coherence = coh,
            direction = dirn,
            driftRate = shiftCols / (2 * h) * tInc / wInc,
            gradEnergy = NA_real_,
            peakCorrelation = peakR,
            peakLag = peakLag,
            peakSharpness = peakSharp
        )
        attr(out, "wccSymmetryParams") <- params
        return(out)
    }

    # Gradients in index units: central differences, one-sided at edges.
    diffAlong <- function(x) {
        n <- length(x)
        d <- numeric(n)
        d[2:(n - 1)] <- (x[3:n] - x[1:(n - 2)]) / 2
        d[1] <- x[2] - x[1]
        d[n] <- x[n] - x[n - 1]
        d
    }
    gT <- apply(m, 2, diffAlong)        # dr/d(row)
    gL <- t(apply(m, 1, diffAlong))     # dr/d(lag column)

    # Zero-change direction vector in (t, lag) is (gL, -gT).
    # Doubled-angle components weighted by |g|^2:
    #   |g|^2 cos(2 phi) = gL^2 - gT^2,  |g|^2 sin(2 phi) = -2 gL gT
    mag2 <- gT^2 + gL^2
    if (weight == "gradient") {
        w <- matrix(1, nT, nL)
        cA <- gL^2 - gT^2
        sA <- -2 * gL * gT
        den <- mag2
    } else {
        ok <- sqrt(mag2) >= minGrad
        w <- ok * 1
        mg2 <- ifelse(ok, mag2, 1)
        cA <- (gL^2 - gT^2) / mg2
        sA <- (-2 * gL * gT) / mg2
        den <- w
    }
    bad <- !is.finite(cA) | !is.finite(sA) | !is.finite(den)
    cA[bad] <- 0; sA[bad] <- 0; den[bad] <- 0; w[bad] <- 0

    C <- rowSums(w * cA)
    S <- rowSums(w * sA)
    D <- rowSums(w * den)
    R <- ifelse(D > 0, sqrt(C^2 + S^2) / D, NA_real_)
    dirn <- ifelse(D > 0, 0.5 * atan2(S, C), NA_real_)

    mag2[bad] <- NA
    out <- data.frame(
        window = rowIdx,
        elapsedTime = (offset + (rowIdx - 1) * wInc) / samplespersecond,
        dispersion = 1 - R,
        coherence = R,
        direction = dirn,
        driftRate = tan(dirn) * tInc / wInc,
        gradEnergy = rowMeans(mag2, na.rm=TRUE),
        peakCorrelation = peakR,
        peakLag = peakLag,
        peakSharpness = peakSharp
    )
    attr(out, "wccSymmetryParams") <- params
    out
}

# ----------------------------------
# Internal: data-driven threshold.  Splits log dispersion into two groups
# (symmetric and not) by maximizing the between-group variance (Otsu's
# method).  Log scale is used because dispersion in symmetric windows is
# often orders of magnitude smaller than in non-symmetric windows, and the
# absolute levels depend on the noise in the data.

iwccOtsuThreshold <- function(d) {
    v <- log(pmax(d[is.finite(d)], 1e-6))
    if (length(v) < 4 || diff(range(v)) == 0) return(exp(max(v, log(1e-6))))
    cuts <- sort(unique(stats::quantile(v, seq(.02, .98, by=.01), names=FALSE)))
    best <- cuts[1]
    bestScore <- -Inf
    for (cut in cuts) {
        a <- v[v <= cut]
        b <- v[v > cut]
        if (length(a) < 2 || length(b) < 2) next
        score <- length(a) * length(b) * (mean(a) - mean(b))^2
        if (score > bestScore) {
            bestScore <- score
            best <- cut
        }
    }
    exp(best)
}

# ----------------------------------
# Internal: delays, in elapsed-time units, between a change in the time
# series and its detection in the windows.  Calibrated on simulated
# dyads with known epochs (see wccSymmetryEpochs).

iwccEpochOffsets <- function(p) {
    if (is.null(p) || is.na(p$wMax) || is.na(p$tMax)) return(c(start=0, end=0))
    sps <- p$samplespersecond
    if (p$method == "shift") {
        h <- p$timeScale * p$wInc
        c(start=(3.3 + 0.70 * p$wMax + 0.78 * h) / sps,
          end=(5.6 + 0.47 * p$wMax - 0.78 * h) / sps)
    } else {
        c(start=0, end=0)
    }
}

# ----------------------------------
# Symmetry formation and breaking epochs.

wccSymmetryEpochs <- function(x, threshold=NA, minDuration=NA,
                              smoothWidth=NA, maxGap=NA,
                              driftTolerance=NULL, ...) {
    if (!is.data.frame(x)) x <- wccZeroChangeDirections(x, ...)
    n <- nrow(x)

    # Defaults tuned on simulated dyads with known epochs; for
    # method="shift" they scale with the comparison distance h (windows).
    p <- attr(x, "wccSymmetryParams")
    isShift <- !is.null(p) && p$method == "shift"
    h <- if (isShift) p$timeScale else NA
    if (identical(threshold, NA)) threshold <- if (isShift) 0.6 else "auto"
    if (is.na(smoothWidth)) smoothWidth <- if (isShift) h + (h %% 2 == 0) else 5
    if (is.na(minDuration)) minDuration <- if (isShift) round(4 * h / 3) else 5
    if (is.na(maxGap)) maxGap <- if (isShift) round(2 * h / 3) else 0

    d <- x$dispersion
    if (smoothWidth > 1 && n >= 3) {
        k <- min(smoothWidth + (smoothWidth %% 2 == 0), n - (n %% 2 == 0))
        dd <- d
        dd[is.na(dd)] <- 1
        d <- as.numeric(stats::runmed(dd, k, endrule="median"))
    }
    if (identical(threshold, "auto")) {
        threshold <- iwccOtsuThreshold(d)
    }
    if (!is.numeric(threshold) || length(threshold) != 1 || !is.finite(threshold)) {
        stop("Warning: threshold must be \"auto\" or a single number between 0 and 1.")
    }
    sym <- !is.na(d) & d <= threshold

    # Fill short interior gaps between symmetric runs.
    if (maxGap > 0) {
        r <- rle(sym)
        ends <- cumsum(r$lengths)
        for (i in seq_along(r$values)) {
            if (!r$values[i] && r$lengths[i] <= maxGap &&
                i > 1 && i < length(r$values)) {
                sym[(ends[i] - r$lengths[i] + 1):ends[i]] <- TRUE
            }
        }
    }

    # Drop runs shorter than minDuration.
    r <- rle(sym)
    r$values[r$values & r$lengths < minDuration] <- FALSE
    sym <- inverse.rle(r)

    r <- rle(sym)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1
    idx <- which(r$values)
    starts <- starts[idx]
    ends <- ends[idx]

    # Optionally split a run where its drift rate changes, so a stationary
    # epoch followed directly by a drifting one is two epochs: one
    # symmetry breaks as another forms.
    drift <- x$driftRate
    if (!is.null(driftTolerance) && length(starts)) {
        dk <- max(3, smoothWidth + (smoothWidth %% 2 == 0))
        newS <- integer(0); newE <- integer(0)
        for (j in seq_along(starts)) {
            s <- starts[j]; e <- ends[j]
            seg <- drift[s:e]
            seg[!is.finite(seg)] <- stats::median(seg[is.finite(seg)])
            if (length(seg) >= dk) seg <- as.numeric(stats::runmed(seg, dk))
            cur <- s; ref <- seg[1]
            for (i in seq_along(seg)[-1]) {
                pos <- s + i - 1
                if (abs(seg[i] - ref) > driftTolerance &&
                    pos - cur >= minDuration && e - pos + 1 >= minDuration) {
                    newS <- c(newS, cur); newE <- c(newE, pos - 1)
                    cur <- pos
                }
                ref <- stats::median(seg[(cur - s + 1):i])
            }
            newS <- c(newS, cur); newE <- c(newE, e)
        }
        starts <- newS; ends <- newE
    }

    # An epoch that runs into the unscored windows at the end of the
    # series (the last timeScale windows of method="shift") has not been
    # seen to break: it is open-ended.
    lastScored <- if (any(is.finite(x$dispersion))) max(which(is.finite(x$dispersion))) else n
    open <- ends >= lastScored - 1

    epochs <- data.frame(
        formationWindow = as.integer(starts),
        breakingWindow = ifelse(ends < n & !open, as.integer(ends + 1), NA_integer_)
    )
    epochs$formationTime <- x$elapsedTime[epochs$formationWindow]
    epochs$breakingTime <- ifelse(is.na(epochs$breakingWindow), NA_real_,
                                  x$elapsedTime[pmin(epochs$breakingWindow, n)])
    epochs$duration <- as.integer(ends - starts + 1)

    # Estimated start and end of the epoch in the time series itself.
    # A window's elapsed time is the end of its data, so symmetry is first
    # seen some time after it forms and is lost soon after it breaks.
    off <- iwccEpochOffsets(attr(x, "wccSymmetryParams"))
    lastTime <- x$elapsedTime[n]
    epochs$startTime <- pmax(epochs$formationTime - off["start"], 0)
    epochs$endTime <- ifelse(is.na(epochs$breakingTime), lastTime,
                             pmin(epochs$breakingTime - off["end"], lastTime))
    # Adjacent epochs separated by a short transition can overlap once the
    # delays are removed; split the overlap at its midpoint.
    sampleStep <- if (is.null(attr(x, "wccSymmetryParams"))) 1 else
        1 / attr(x, "wccSymmetryParams")$samplespersecond
    for (i in seq_len(nrow(epochs))[-1]) {
        if (epochs$startTime[i] <= epochs$endTime[i - 1]) {
            mid <- (epochs$startTime[i] + epochs$endTime[i - 1]) / 2
            epochs$endTime[i - 1] <- mid
            epochs$startTime[i] <- mid + sampleStep
        }
    }

    seg <- function(v, f) vapply(seq_along(starts), function(i)
        f(v[starts[i]:ends[i]]), 0)
    epochs$meanDispersion <- seg(x$dispersion, function(v) mean(v, na.rm=TRUE))
    epochs$medianDriftRate <- seg(drift, function(v) stats::median(v, na.rm=TRUE))
    if (!is.null(x$peakCorrelation)) {
        epochs$peakCorrelation <- seg(x$peakCorrelation, function(v) stats::median(v, na.rm=TRUE))
        epochs$peakLag <- seg(x$peakLag, function(v) stats::median(v, na.rm=TRUE))
        lagStep <- if (is.null(attr(x, "wccSymmetryParams"))) 1 else
            with(attr(x, "wccSymmetryParams"), tInc / samplespersecond)
        epochs$peakLagStability <- seg(x$peakLag, function(v)
            mean(abs(v - stats::median(v, na.rm=TRUE)) <= lagStep + 1e-9, na.rm=TRUE))
        epochs$peakSharpness <- seg(x$peakSharpness, function(v) stats::median(v, na.rm=TRUE))
    }

    x$smoothedDispersion <- d
    x$symmetric <- sym
    list(epochs=epochs, windows=x, threshold=threshold)
}
