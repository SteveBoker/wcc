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
# Program: wccPlot.R
#  Author: Steven Boker
#    Date: Wed Jan 28 09:05:29 EST 2026
#
# This function takes the filename of a wcc output csv file created by wccAggregate() and
#   plots a wcc graphic.
#
# ---------------------------------------------------------------------
# Revision History
#  Steve Boker  -- Wed Jan 28 09:05:30 EST 2026
#      Created wccPlot.R
# ---------------------------------------------------------------------
# Revision History
#  Steve Boker  -- Wed Apr 8 09:44:06 EDT 2026
#      Changed function so that wccPlot calls wccCalc and wccPeakPick rather than reading
#      this information from disk.
#      Also removed the pdffile option.  
#
# ---------------------------------------------------------------------

library(pheatmap)
library(gtable)

wccPlot <- function(inSeries1=NA, inSeries2=NA, startwindow=1, endwindow=200, wMax=50, tMax=50, wInc=1, tInc=1, Lsize=8, pspan=.25, type="Max", samplespersecond=1, method=c("c", "cumr", "cumc", "r"), ...) {

    # Deprecation: allow old windcross argument
    dots <- list(...)
    if ("windcross" %in% names(dots)) {
        warning("Argument 'windcross' is deprecated; use method=\"c\" or method=\"r\" instead.")
        method <- if (isTRUE(dots$windcross)) "c" else "r"
    } else {
        method <- match.arg(method)
    }

    if (!is.numeric(inSeries1) | !is.numeric(inSeries2) | !is.vector(inSeries1) | !is.vector(inSeries2) | length(inSeries1) != length(inSeries2)) {
        stop(paste0("Warning: inSeries1 and inSeries2 must be numeric vectors of equal length."))
    }
    if (!is.numeric(wMax) | wMax < 5  ) {
        stop(paste0("Warning: wMax must be a numeric greater than 5."))
    }
    if ( !is.numeric(wInc) |  wInc < 1 ) {
        stop(paste0("Warning: wInc must be a numeric greater than or equal to 1."))
    }
    #    tFrame <- read.csv(savefile)
    
    # ----------------------------------
    # calculate the windowed cross correlation of inSeries1 and inSeries2.

    wccMatrix <- wccCalc(inSeries1, inSeries2, wMax=wMax, tMax=tMax, wInc=wInc, tInc=tInc, method=method)
    wccZeroRow <- floor((dim(wccMatrix)[2] + 1) / 2)
    wccColumns <- dim(wccMatrix)[1]
    wccRows <- dim(wccMatrix)[2]

    # ----------------------------------
    # Check whether startwindow and endwindow are possible given wccColumns.

    if (startwindow < 1 ) {
        stop(paste0("Warning: startwindow must be greater than or equal to 1"))
    }
    if (startwindow < 1 | startwindow > wccColumns ) {
        stop(paste0("Warning: startwindow must be less than the total number of WCC windows, ", wccColumns))
    }

    if (endwindow < startwindow ) {
        stop(paste0("Warning: endwindow must be greater than startwindow"))
    }
    if (startwindow < 1 | startwindow > wccColumns ) {
        stop(paste0("Warning: endwindow must be less than or equal to the total number of WCC windows, ", wccColumns))
    }

    ppOut <- wccPeakPick(wccMatrix, Lsize=Lsize, pspan=pspan, type=type)
    ppOut$maxIndex <- ppOut$maxIndex*.5

    # ----------------------------------
    # Extract the matrix of correlations into the format desired by image().

    plotColumns <- (endwindow - startwindow) + 1
    theX <- (wMax + tMax + seq(startwindow, endwindow, by=1)) * wInc / samplespersecond
    theY <- seq(-floor(wccRows/2),floor(wccRows/2), by=1) * tInc / samplespersecond
    theWCC <- matrix(NA, length(theY), length(theX))
    i <- 1
    for (j in startwindow:endwindow) {
        for (k in 1:wccRows) {
            theWCC[(wccRows + 1) - k, i] <- wccMatrix[j, k]
        }
        i <- i+1
    }
    # ----------------------------------
    # Plot results.
    
    oldpar <- par(no.readonly = TRUE) # code line i
    on.exit(par(oldpar)) # code line i + 1

    # Layout: main plot + legend
    p <- layout(matrix(c(1,2), ncol = 2), widths = c(6,1))

    par(mar = c(5, 4, 2, 1))
    image(theX, theY, t(theWCC), col=heat.colors(20, alpha=1, rev = FALSE),
         xlab="Elapsed Time (seconds)",
         ylab="Lag Offset (seconds)",
         cex.lab=1.75, cex.axis=1.5, cex.main=1.75, 
         cex=1.75, mgp=c(2.5,.75,0),zlim=c(-1,1))
    lines((wMax + tMax + c(startwindow:endwindow) * wInc) / samplespersecond, 
          -ppOut$maxIndex[startwindow:endwindow] * tInc / samplespersecond, 
          type='l', lwd=2)
    lines(c((wMax + tMax + startwindow * wInc) / samplespersecond, (wMax + tMax + endwindow * wInc) / samplespersecond), c(0,0), type='l', lty=2, lwd=2)

    par(mar = c(5, 1, 2, 3))  # adjust margins
    image(x = 1,
          y = seq(-1,1, length.out=20),
          z = matrix(seq(-1,1, length.out=20), nrow= 1),
          xlab=" ",
          ylab="",
          col = heat.colors(20, alpha=1, rev = FALSE),
          axes = FALSE)

    axis(4)  # add axis on right
}