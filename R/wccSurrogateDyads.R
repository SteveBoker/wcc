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
# Program: wccSurrogateDyads.R
#  Author: Steven Boker
#    Date: Wed Jan 28 09:08:34 EST 2026
#
# This function takes an array of timeseries, shuffles them, runs wccAggregate()
#   on each shuffled pair and returns an array of wcc objects that represent a
#   distribution of the null hypothesis that pairs of timeseries are unrelated.
#
# ---------------------------------------------------------------------
# Revision History
#  Steve Boker  -- Wed Jan 28 09:08:36 EST 2026
#      Created wccSurrogateDyads.R
#
# ---------------------------------------------------------------------




wccSurrogateDyads <- function(inArray1=NA, inArray2=NA, wMax=50, tMax=50, wInc=1, tInc=1, Lsize=8, pspan=.25, type="Max", nSurrogates=NA, method=c("c", "cumr", "cumc", "cumcuda", "r"), embedD=9, ...) {
    # Deprecation: allow old windcross argument
    dots <- list(...)
    if ("windcross" %in% names(dots)) {
        warning("Argument 'windcross' is deprecated; use method=\"c\" or method=\"r\" instead.")
        method <- if (isTRUE(dots$windcross)) "c" else "r"
    } else {
        method <- match.arg(method)
    }
    #Note that wccAggNames must match the definition in wccAggregate, wccFindDyadParam and wccSurrogateDyads
    wccAggNames <- c("wMax", "tMax", "wInc", "tInc", "Lsize", "pspan", "type", "samples", "windows", 
                    "pctMissing", "pctMissingWindows", "maxMean", "maxVar", "totalMean", "totalVar", "zeroLagMean", 'zeroLagVar', 
                    "lagMean", "lagVar", "dlagMean", "dlagVar", "d2lagMean", "d2lagVar", "elapsedSeconds")
    if (!is.numeric(inArray1) | !is.numeric(inArray2) | !is.matrix(inArray1) | !is.matrix(inArray2) | dim(inArray1)[1] != dim(inArray2)[1] | dim(inArray1)[2] != dim(inArray2)[2]) {
        stop(paste0("Warning: inArray1 and inArray2 must be numeric matrices of the same dimension."))
    }
    if (!is.numeric(embedD) | embedD < 5  ) {
        stop(paste0("Warning: embedD must be a numeric greater than or equal to 5."))
    }
    if (!is.numeric(wMax) | wMax < 5  ) {
        stop(paste0("Warning: wMax must be a numeric greater than 5."))
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
    surrogateFrame <- data.frame(matrix(NA, nrow=nSurrogates, ncol=length(wccAggNames))) 
    names(surrogateFrame) <- wccAggNames
    for (i in 1:nSurrogates) {
        j <- floor(runif(1, min=1, max=dim(inArray1)[1]))
        k <- j
        while(k == j) {
            k <- floor(runif(1, min=1, max=dim(inArray1)[1]))            
        }
        message(paste0("Iteration ", i, " of ", nSurrogates))
        surrogateFrame[i,] <- wccAggregate(inSeries1=inArray1[j,], inSeries2=inArray2[k,], wMax=wMax, tMax=tMax, wInc=wInc, 
            tInc=tInc, Lsize=Lsize, pspan=pspan, type=type, method=method, embedD=embedD)
    }
    return(surrogateFrame)
}
