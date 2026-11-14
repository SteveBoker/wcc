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
# Program: wccFindDyadParam.R
#  Author: Steven Boker
#    Date: Mon Feb 2 13:30:58 EST 2026
#
# This function receives a set of parameters to optimize and two arrays of
#   timeseries.  It repeatedly calls wccSurrogateDyads() and wccAggregate() 
#   for each parameter choice and returns the parameters for which real results are
#   most different than surrogate distributions for the results defined by opt1, opt2, and opt3.
#
# ---------------------------------------------------------------------
# Revision History
#  Steve Boker  -- Mon Feb 2 13:31:03 EST 2026
#      Created wccFindDyadParam.R
#
# ---------------------------------------------------------------------


wccFindDyadParam <- function(inArray1=NA, inArray2=NA, wMaxvector=c(50), tMaxvector=c(50), wIncvector=c(1), tIncvector=c(1),
                             Lsizevector=c(8), pspanvector=c(.25), type="Max", nSurrogates=NA, samplespersecond=1, method=c("c", "cumr", "cumc", "cumcuda", "r"), embedD=9, ...) {
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
    if (!is.numeric(wMaxvector) | min(wMaxvector) < 5  ) {
        stop(paste0("Warning: wMax must be a numeric greater than 5."))
    }
    if ( !is.numeric(wIncvector) |  min(wIncvector) < 1 ) {
        stop(paste0("Warning: wInc must be a numeric greater than or equal to 1."))
    }
    if (!is.numeric(tIncvector) | min(tIncvector) < 1  ) {
        stop(paste0("Warning: tInc must be a numeric greater than or equal to 1."))
    }
    if (!is.numeric(tMaxvector) | min(tMaxvector) < 5  ) {
        stop(paste0("Warning: tMax must be a numeric greater than 5."))
    }
    if (!is.numeric(Lsizevector) | min(Lsizevector) < 5  ) {
        stop(paste0("Warning: Lsize must be a numeric greater than 5."))
    }
    if (!is.numeric(pspanvector) | min(pspanvector) < .05  ) {
        stop(paste0("Warning: pspan must be a numeric greater than 0.05."))
    }
    wccTestNames <- c(wccAggNames, "maxMeanKS", "maxVarKS", "totalMeanKS", "totalVarKS", "zeroLagMeanKS", 'zeroLagVarKS', "lagMeanKS", "lagVaKSr", "dlagMeanKS", "dlagVaKSr", "d2lagMeanKS", "d2lagVaKSr",
     "maxMeanQdiff", "maxVarQdiff", "totalMeanQdiff", "totalVarQdiff", "zeroLagMeanQdiff", 'zeroLagVarQdiff', "lagMeanQdiff", "lagVarQdiff", "dlagMeanQdiff", "dlagVarQdiff", "d2lagMeanQdiff", "d2lagVarQdiff")

    totalTests <- length(wIncvector) * length(wMaxvector) * length(tMaxvector) * length(tIncvector) * length(Lsizevector) * length(pspanvector)
    testFrame <- data.frame(matrix(NA, nrow=totalTests, ncol=length(wccAggNames)))
    names(testFrame) <- wccAggNames
    nDyads <- dim(inArray1)[1]

    # Batched path: surrogate pairings are drawn once and reused for every
    # parameter combination; real dyads and surrogates are computed together
    # in a single wccCalcBatch call per (wInc, wMax, tMax, tInc) combination.
    useBatch <- method %in% c("cumc", "cumcuda")
    if (useBatch) {
        if (anyNA(inArray1) || anyNA(inArray2)) {
            stop(paste0("Warning: method \"", method, "\" does not support missing data. Use method=\"c\" or method=\"r\"."))
        }
        surrogatePairs <- matrix(NA_integer_, nrow=nSurrogates, ncol=2)
        for (s in 1:nSurrogates) {
            j <- sample.int(nDyads, 1)
            k <- j
            while (k == j) {
                k <- sample.int(nDyads, 1)
            }
            surrogatePairs[s,] <- c(j, k)
        }
        allPairs <- rbind(cbind(1:nDyads, 1:nDyads), surrogatePairs)
    }

    testIndex <- 1
    for(testwInc in wIncvector) {
        for(testwMax in wMaxvector) {
            for(testtMax in tMaxvector) {
                for(testtInc in tIncvector) {
                    if (useBatch) {
                        batchStart <- Sys.time()
                        wccGrids <- wccCalcBatch(inArray1, inArray2, pairs=allPairs, wMax=testwMax, tMax=testtMax,
                                wInc=testwInc, tInc=testtInc, method=method)
                    }
                    for(testLsize in Lsizevector) {
                        for(testpspan in pspanvector) {
                            message(paste0("Parameter Test ", testIndex, " of ", totalTests))
                            realFrame <- data.frame(matrix(NA, nrow=dim(inArray1)[1], ncol=length(wccAggNames)))
                            names(realFrame) <- wccAggNames
                            if (useBatch) {
                                aggFrame <- do.call(rbind, lapply(seq_len(dim(wccGrids)[3]), function(p) {
                                    wccAggregateGrid(wccGrids[,,p], wMax=testwMax, tMax=testtMax, wInc=testwInc, tInc=testtInc,
                                        Lsize=testLsize, pspan=testpspan, type=type, samplespersecond=samplespersecond,
                                        embedD=embedD, nSamples=dim(inArray1)[2], pctmissing=0, startTime=batchStart)
                                }))
                                names(aggFrame) <- wccAggNames
                                realFrame[,] <- aggFrame[1:nDyads,]
                                surrogateFrame <- aggFrame[(nDyads+1):(nDyads+nSurrogates),]
                            }
                            else {
                                surrogateFrame <- wccSurrogateDyads(inArray1=inArray1, inArray2=inArray2, wInc=testwInc, wMax=testwMax, tMax=testtMax,
                                        tInc=testtInc, Lsize=testLsize, pspan=testpspan, type=type, nSurrogates=nSurrogates, method=method, embedD=embedD)
                            }
                            for(i in 1:dim(inArray1)[1]) {
                                if (!useBatch) {
                                    realFrame[i,] <- wccAggregate(inSeries1=inArray1[i,], inSeries2=inArray2[i,], wInc=testwInc, wMax=testwMax, tMax=testtMax,
                                        tInc=testtInc, Lsize=testLsize, pspan=testpspan, type=type, method=method)
                                }
                                realFrame$maxMean[i] <- sum(surrogateFrame$maxMean > realFrame$maxMean[i]) / nSurrogates
                                realFrame$maxVar[i] <- sum(surrogateFrame$maxVar > realFrame$maxVar[i]) / nSurrogates
                                realFrame$totalMean[i] <- sum(surrogateFrame$totalMean > realFrame$totalMean[i]) / nSurrogates
                                realFrame$totalVar[i] <- sum(surrogateFrame$totalVar > realFrame$totalVar[i]) / nSurrogates
                                realFrame$zeroLagMean[i] <- sum(surrogateFrame$zeroLagMean > realFrame$zeroLagMean[i]) / nSurrogates
                                realFrame$zeroLagVar[i] <- sum(surrogateFrame$zeroLagVar > realFrame$zeroLagVar[i]) / nSurrogates
                                realFrame$lagMean[i] <- sum(surrogateFrame$lagMean > realFrame$lagMean[i]) / nSurrogates
                                realFrame$lagVar[i] <- sum(surrogateFrame$lagVar > realFrame$lagVar[i]) / nSurrogates
                                realFrame$dlagMean[i] <- sum(surrogateFrame$dlagMean > realFrame$lagMean[i]) / nSurrogates
                                realFrame$dlagVar[i] <- sum(surrogateFrame$dlagVar > realFrame$lagVar[i]) / nSurrogates
                            }
                            testFrame$wInc[testIndex] <- testwInc
                            testFrame$wMax[testIndex] <- testwMax
                            testFrame$tMax[testIndex] <- testtMax
                            testFrame$tInc[testIndex] <- testtInc
                            testFrame$Lsize[testIndex] <- testLsize
                            testFrame$pspan[testIndex] <- testpspan
                            tPPtype <- 1
                            if (type == "Min") tPPtype <- -1
                            testFrame$type[testIndex] <- tPPtype
                            testFrame$samples[testIndex] <- realFrame$samples[1]
                            testFrame$windows[testIndex] <- realFrame$windows[1]
                            testFrame$pctMissing[testIndex] <- mean(realFrame$pctMissing)
                            testFrame$pctMissingWindows[testIndex] <- mean(realFrame$pctMissingWindows)
                            testFrame$maxMean[testIndex] <- sum(realFrame$maxMean>.95 | realFrame$maxMean<.05) / dim(inArray1)[1]
                            testFrame$maxVar[testIndex] <- sum(realFrame$maxVar>.95 | realFrame$maxVar<.05) / dim(inArray1)[1]
                            testFrame$totalMean[testIndex] <- sum(realFrame$totalMean>.95 | realFrame$totalMean<.05) / dim(inArray1)[1]
                            testFrame$totalVar[testIndex] <- sum(realFrame$totalVar>.95 | realFrame$totalVar<.05) / dim(inArray1)[1]
                            testFrame$zeroLagMean[testIndex] <- sum(realFrame$zeroLagMean>.95 | realFrame$zeroLagMean<.05) / dim(inArray1)[1]
                            testFrame$zeroLagVar[testIndex] <- sum(realFrame$zeroLagVar>.95 | realFrame$zeroLagVar<.05) / dim(inArray1)[1]
                            testFrame$lagMean[testIndex] <- sum(realFrame$lagMean>.95 | realFrame$lagMean<.05) / dim(inArray1)[1]
                            testFrame$lagVar[testIndex] <- sum(realFrame$lagVar>.95 | realFrame$lagVar<.05) / dim(inArray1)[1]
                            testFrame$dlagMean[testIndex] <- sum(realFrame$dlagMean>.95 | realFrame$dlagMean<.05) / dim(inArray1)[1]
                            testFrame$dlagVar[testIndex] <- sum(realFrame$dlagVar>.95 | realFrame$dlagVar<.05) / dim(inArray1)[1]
                            testFrame$d2lagMean[testIndex] <- sum(realFrame$d2lagMean>.95 | realFrame$d2lagMean<.05) / dim(inArray1)[1]
                            testFrame$d2lagVar[testIndex] <- sum(realFrame$d2lagVar>.95 | realFrame$d2lagVar<.05) / dim(inArray1)[1]
                            testFrame$elapsedSeconds[testIndex] <- sum(surrogateFrame$elapsedSeconds) + sum(realFrame$elapsedSeconds)
                            testFrame$maxMeanKS[testIndex] <- ks.test(realFrame$maxMean, surrogateFrame$maxMean, alternative="two.sided")$p.value
                            testFrame$maxVarKS[testIndex] <- ks.test(realFrame$maxVar, surrogateFrame$maxVar, alternative="two.sided")$p.value
                            testFrame$totalMeanKS[testIndex] <- ks.test(realFrame$totalMean, surrogateFrame$totalMean, alternative="two.sided")$p.value
                            testFrame$totalVarKS[testIndex] <- ks.test(realFrame$totalVar, surrogateFrame$totalVar, alternative="two.sided")$p.value
                            testFrame$zeroLagMeanKS[testIndex] <- ks.test(realFrame$zeroLagMean, surrogateFrame$zeroLagMean, alternative="two.sided")$p.value
                            testFrame$zeroLagVarKS[testIndex] <- ks.test(realFrame$zeroLagVar, surrogateFrame$zeroLagVar, alternative="two.sided")$p.value
                            testFrame$lagMeanKS[testIndex] <- ks.test(realFrame$lagMean, surrogateFrame$lagMean, alternative="two.sided")$p.value
                            testFrame$lagVarKS[testIndex] <- ks.test(realFrame$lagVar, surrogateFrame$lagVar, alternative="two.sided")$p.value
                            testFrame$dlagMeanKS[testIndex] <- ks.test(realFrame$dlagMean, surrogateFrame$dlagMean, alternative="two.sided")$p.value
                            testFrame$dlagVarKS[testIndex] <- ks.test(realFrame$dlagVar, surrogateFrame$dlagVar, alternative="two.sided")$p.value
                            testFrame$d2lagMeanKS[testIndex] <- ks.test(realFrame$d2lagMean, surrogateFrame$d2lagMean, alternative="two.sided")$p.value
                            testFrame$d2lagVarKS[testIndex] <- ks.test(realFrame$d2lagVar, surrogateFrame$d2lagVar, alternative="two.sided")$p.value

                            testFrame$maxMeanQdiff[testIndex] <- mean(quantile(realFrame$maxMean, probs=seq(.1,.9,by=.1), na.rm=TRUE)- quantile(surrogateFrame$maxMean, probs=seq(.1,.9,by=.1), na.rm=TRUE))
                            testFrame$maxVarQdiff[testIndex] <- mean(quantile(realFrame$maxVar, probs=seq(.1,.9,by=.1), na.rm=TRUE)- quantile(surrogateFrame$maxVar, probs=seq(.1,.9,by=.1), na.rm=TRUE))
                            testFrame$totalMeanQdiff[testIndex] <- mean(quantile(realFrame$totalMean, probs=seq(.1,.9,by=.1), na.rm=TRUE)- quantile(surrogateFrame$totalMean, probs=seq(.1,.9,by=.1), na.rm=TRUE))
                            testFrame$totalVarQdiff[testIndex] <- mean(quantile(realFrame$totalVar, probs=seq(.1,.9,by=.1), na.rm=TRUE)- quantile(surrogateFrame$totalVar, probs=seq(.1,.9,by=.1), na.rm=TRUE))
                            testFrame$zeroLagMeanQdiff[testIndex] <- mean(quantile(realFrame$zeroLagMean, probs=seq(.1,.9,by=.1), na.rm=TRUE)- quantile(surrogateFrame$zeroLagMean, probs=seq(.1,.9,by=.1), na.rm=TRUE))
                            testFrame$zeroLagVarQdiff[testIndex] <- mean(quantile(realFrame$zeroLagVar, probs=seq(.1,.9,by=.1), na.rm=TRUE)- quantile(surrogateFrame$zeroLagVar, probs=seq(.1,.9,by=.1), na.rm=TRUE))
                            testFrame$lagMeanQdiff[testIndex] <- mean(quantile(realFrame$lagMean, probs=seq(.1,.9,by=.1), na.rm=TRUE)- quantile(surrogateFrame$lagMean, probs=seq(.1,.9,by=.1), na.rm=TRUE))
                            testFrame$lagVarQdiff[testIndex] <- mean(quantile(realFrame$lagVar, probs=seq(.1,.9,by=.1), na.rm=TRUE)- quantile(surrogateFrame$lagVar, probs=seq(.1,.9,by=.1), na.rm=TRUE))
                            testFrame$dlagMeanQdiff[testIndex] <- mean(quantile(realFrame$dlagMean, probs=seq(.1,.9,by=.1), na.rm=TRUE)- quantile(surrogateFrame$dlagMean, probs=seq(.1,.9,by=.1), na.rm=TRUE))
                            testFrame$dlagVarQdiff[testIndex] <- mean(quantile(realFrame$dlagVar, probs=seq(.1,.9,by=.1), na.rm=TRUE)- quantile(surrogateFrame$dlagVar, probs=seq(.1,.9,by=.1), na.rm=TRUE))
                            testFrame$d2lagMeanQdiff[testIndex] <- mean(quantile(realFrame$d2lagMean, probs=seq(.1,.9,by=.1), na.rm=TRUE)- quantile(surrogateFrame$d2lagMean, probs=seq(.1,.9,by=.1), na.rm=TRUE))
                            testFrame$d2lagVarQdiff[testIndex] <- mean(quantile(realFrame$d2lagVar, probs=seq(.1,.9,by=.1), na.rm=TRUE)- quantile(surrogateFrame$d2lagVar, probs=seq(.1,.9,by=.1), na.rm=TRUE))
                            testIndex <- testIndex + 1 
                        }        
                    }    
                }  
            }
        } 
    }
    return(testFrame)
}