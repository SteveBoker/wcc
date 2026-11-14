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
# Program: wccAggregate.R
#  Author: Steven Boker
#    Date: Wed Jan 28 08:55:03 EST 2026
#
# This function takes a pair of timeseries runs wcc and peak picking and returns 
#   a wccAgg dataframe with aggregated statistics.  Optionally it saves the wccCalc matrix
#   and peakpicking vector to a csv file.
#
#
# ---------------------------------------------------------------------
# Revision History
#  Steve Boker  -- Wed Jan 28 08:55:06 EST 2026
#      Created wccAggregate.R
#
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# wccAggregate returns a data frame with one row and the following columns
#     wInc            --  Window increment in samples
#     wMax            --  Window size in samples
#     tMax            --  Maximum lag in samples
#     tInc            --  Lag increment in samples
#     Lsize            --  Size of bounding box for finding peaks/valleys
#     pspan            --  Smoothing parameter from LOESS
#     type             --  "Max" = find peaks, "Min"= find valleys
#     samples            --  Total samples in time series
#     windows            --  Total windows after running WCC
#     pctMissing         --  Proportion of missing values in the two time series
#     pctMissingWindows  --  Proportion of missing windows after running WCC
#     maxMean            --  Mean of the maximum WCC closest to zero lag
#     maxVar             --  Variance of the maximum WCC closest to zero lag
#     totalMean          --  Mean of all values in the WCC matrix
#     totalVar           --  Variance of all values in the WCC matrix
#     zeroLagMean        --  Mean of the values in the WCC matrix at zero lag
#     zeroLagVar         --  Variance of the values in the WCC matrix at zero lag
#     lagMean            --  Mean of the lag found by peak picking
#     lagVar             --  Variance of the lag found by peak picking
#     dlagMean           --  Mean of the first derivative of the lag found by peak picking
#     dlagVar            --  Variance of the first derivative of the lag found by peak picking
#     elapsedTime        --  Elapsed time in seconds
# ---------------------------------------------------------------------


wccAggregate <- function(inSeries1=NA, inSeries2=NA, wMax=50, tMax=50, wInc=1, tInc=1, Lsize=8, pspan=.25, type="Max", samplespersecond=1, method=c("c", "cumr", "cumc", "cumcuda", "r"), embedD=9, ...) {

    # Deprecation: allow old windcross argument
    dots <- list(...)
    if ("windcross" %in% names(dots)) {
        warning("Argument 'windcross' is deprecated; use method=\"c\" or method=\"r\" instead.")
        method <- if (isTRUE(dots$windcross)) "c" else "r"
    } else {
        method <- match.arg(method)
    }

    #Note that wccAggNames must match the definition in wccFindDyadParam and wccSurrogateDyads
    wccAggNames <- c("wMax", "tMax", "wInc", "tInc", "Lsize", "pspan", "type", "samples", "windows", 
                    "pctMissing", "pctMissingWindows", "maxMean", "maxVar", "totalMean", "totalVar", "zeroLagMean", 'zeroLagVar', 
                    "lagMean", "lagVar", "dlagMean", "dlagVar", "d2lagMean", "d2lagVar", "elapsedSeconds")

    startTime <- Sys.time()

    if (!is.numeric(inSeries1) | !is.numeric(inSeries2) | !is.vector(inSeries1) | !is.vector(inSeries2) | length(inSeries1) != length(inSeries2)) {
        stop(paste0("Warning: inSeries1 and inSeries2 must be numeric vectors."))
    }
    if (length(inSeries1) != length(inSeries2)) {
        stop(paste0("Warning: inSeries1 and inSeries2 must be vectors of equal length."))
    }
    if (!is.numeric(embedD) | embedD < 5  ) {
        stop(paste0("Warning: embedD must be a numeric greater than or equal to 5."))
    }
    if (!is.numeric(wMax) | wMax < 5  ) {
        stop(paste0("Warning: wMax must be a numeric greater than or equal to 5."))
    }
    if ( !is.numeric(wInc) |  wInc < 1 ) {
        stop(paste0("Warning: wInc must be a numeric greater than or equal to 1."))
    }
    if (!is.numeric(tInc) | tInc < 1  ) {
        stop(paste0("Warning: tInc must be a numeric greater than or equal to 1."))
    }
    if (!is.numeric(tMax) | tMax < 5  ) {
        stop(paste0("Warning: tMax must be a numeric greater than 5."))
    }
    if (!is.numeric(Lsize) | Lsize < 5  ) {
        stop(paste0("Warning: Lsize must be a numeric greater than 5."))
    }
    if (!is.numeric(pspan) | pspan < .05  ) {
        stop(paste0("Warning: tMax must be a numeric greater than 0.05."))
    }
    wccOut <- wccCalc(inSeries1, inSeries2, wMax=wMax, tMax=tMax, wInc=wInc, tInc=tInc, method=method)
    pctmissing <- sum(is.na(c(inSeries1,inSeries2)))/(2*length(inSeries1))
    wccAgg <- wccAggregateGrid(wccOut, wMax=wMax, tMax=tMax, wInc=wInc, tInc=tInc, Lsize=Lsize, pspan=pspan,
            type=type, samplespersecond=samplespersecond, embedD=embedD,
            nSamples=length(inSeries1), pctmissing=pctmissing, startTime=startTime)
    return(wccAgg)
}

# ----------------------------------
# Post-processing of one WCC grid into the aggregate statistics row.
# Shared by wccAggregate (single dyad) and the batched path in
# wccFindDyadParam, which applies it to each slab of a wccCalcBatch result.

wccAggregateGrid <- function(wccOut, wMax, tMax, wInc, tInc, Lsize, pspan, type,
                             samplespersecond, embedD, nSamples, pctmissing,
                             startTime=Sys.time()) {
    pctmissingwindows <- sum(is.na(wccOut[,2]))/dim(wccOut)[1]
    wccOut[is.na(wccOut)] <- 0
    ppOut <- wccPeakPick(wccOut, Lsize=Lsize, pspan=pspan, type=type)
    if (type=="Min") {
        ppOut$maxIndex <- ppOut$minIndex
        ppOut$maxValue <- ppOut$minValue
    }
    ppOut$maxIndex <- ppOut$maxIndex*.5
    PPderivs <- (wccGLLAEmbed(ppOut$maxIndex, embed=embedD, tau=1, idColumn=FALSE) %*% wccGLLAWMatrix(embed=embedD, tau=1, deltaT=samplespersecond, order=2))
    ttype <- 1
    if (type=="Min") {
        ttype <- -1
    }
    wccAgg <- data.frame(
            wInc=wInc,
            wMax=wMax,
            tMax=tMax,
            tInc=tInc,
            Lsize=Lsize,
            pspan=pspan,
            type = ttype,
            samples=nSamples,
            windows=dim(wccOut)[1],
            pctmissing=pctmissing,
            pctmissingwindows=pctmissingwindows,
            maxMean=mean(atanh(c(ppOut$maxValue)), na.rm=TRUE),
            maxVar=var(atanh(c(ppOut$maxValue)), na.rm=TRUE),
            totalMean=mean(abs(atanh(c(wccOut))), na.rm=TRUE),
            totalVar=var(atanh(c(wccOut)), na.rm=TRUE),
            zeroLagMean=mean(atanh(c(wccOut[,tMax+1])), na.rm=TRUE),
            zeroLagVar=var(atanh(c(wccOut[,tMax+1])), na.rm=TRUE),
            lagMean=mean(PPderivs[,1], na.rm=TRUE),
            lagVar=var(PPderivs[,1], na.rm=TRUE),
            dlagMean=mean(PPderivs[,2], na.rm=TRUE),
            dlagVar=var(PPderivs[,2], na.rm=TRUE),
            d2lagMean=mean(PPderivs[,3], na.rm=TRUE),
            d2lagVar=var(PPderivs[,3], na.rm=TRUE),
            elapsedSeconds= Sys.time()-startTime
        )
    return(wccAgg)
}