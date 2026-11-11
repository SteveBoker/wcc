%
%   Copyright 2001-2026 by the individuals mentioned in the source code history
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

\name{wccSymmetryEpochs}
\alias{wccSymmetryEpochs}

\title{wccSymmetryEpochs}

\description{
This function estimates the elapsed times at which symmetry forms and symmetry breaks in a windowed cross correlation matrix.  An epoch of symmetry is a run of consecutive windows in which the dispersion of the directions of zero change, as calculated by \code{\link{wccZeroChangeDirections}}, stays at or below a threshold.  Symmetry formation is the first window of each epoch and symmetry breaking is the first window after it.
}

\usage{
wccSymmetryEpochs(x, threshold=NA, minDuration=NA, smoothWidth=NA,
     maxGap=NA, driftTolerance=NULL, ...)
}

\arguments{
\item{x}{Either the data frame returned by \code{\link{wccZeroChangeDirections}}, or a windowed cross correlation matrix as returned by \code{\link{wccCalc}}.  If a matrix is given, \code{wccZeroChangeDirections} is called on it with the arguments in \code{...}.}
\item{threshold}{NA (default), \sQuote{auto}, or a number.  A window is counted as symmetric when its smoothed dispersion is at or below the threshold.  NA means .6 for \code{method="shift"} (the later pattern correlates at least .4 with the shifted earlier pattern) and \sQuote{auto} for \code{method="gradient"}.  With \sQuote{auto}, the threshold is chosen from the data by splitting the logarithm of the smoothed dispersion into two groups, symmetric and not symmetric, so that the variance between the groups is largest (Otsu's method).  This assumes that the data contain both symmetric and non-symmetric intervals.  \sQuote{auto} is better for clean, nearly deterministic signals; see Details.  The threshold used is returned.}
\item{minDuration}{NA (default) or an integer.  The minimum number of consecutive symmetric windows required for an epoch.  Shorter runs are ignored.  NA means \code{round(4*h/3)} for \code{method="shift"}, where h is the \code{timeScale} of \code{\link{wccZeroChangeDirections}} (20 windows with wMax=30 and wInc=1), and 5 for \code{method="gradient"}.}
\item{smoothWidth}{NA (default) or an integer.  The width, in windows, of a running median applied to the dispersion before it is compared with \code{threshold}.  This suppresses isolated windows that cross the threshold.  NA means h (made odd) for \code{method="shift"} and 5 for \code{method="gradient"}.  Even values are increased by 1.  Set to 1 for no smoothing.}
\item{maxGap}{NA (default) or an integer.  Epochs separated by no more than this many non-symmetric windows are merged into a single epoch.  Gaps are filled before \code{minDuration} is applied.  NA means \code{round(2*h/3)} for \code{method="shift"} and 0 for \code{method="gradient"}.}
\item{driftTolerance}{Either NULL (default) or a numeric.  If NULL, any uniform drift counts as one symmetry, even if the rate of drift changes partway through an epoch.  If a numeric, an epoch is split wherever its drift rate moves more than this amount away from the median drift rate so far in the epoch.  For example, a stationary interval followed directly by an interval of drift is then reported as two epochs, the breaking of one symmetry and the formation of another.  The units are those of \code{driftRate} in \code{\link{wccZeroChangeDirections}}: lag per unit of elapsed time.}
\item{...}{Arguments passed to \code{\link{wccZeroChangeDirections}} when \code{x} is a matrix, for example \code{wMax}, \code{tMax}, \code{wInc}, \code{tInc}, and \code{samplespersecond}.}
}

\details{
The dispersion is first smoothed with a running median of width \code{smoothWidth}.  Windows whose smoothed dispersion is at or below \code{threshold} are marked as symmetric.  Gaps of up to \code{maxGap} windows between symmetric runs are filled, runs shorter than \code{minDuration} windows are removed, and, if \code{driftTolerance} is given, runs are split where their drift rate changes.

Each window of \code{wccCalc} covers \code{wMax} samples ending at its elapsed time, and lagged windows reach back a further \code{tMax} samples.  Symmetry that forms in the time series is therefore first seen in windows some time later, and symmetry that breaks is lost soon after.  \code{startTime} and \code{endTime} correct \code{formationTime} and \code{breakingTime} for these delays.  For \code{method="shift"} the corrections, in samples, are 3.3 + .70 wMax + .78 h wInc for the start and 5.6 + .47 wMax - .78 h wInc for the end.  They were calibrated, so that the median errors are zero, on simulated dyads with known epochs (tMax = 30, wMax from 25 to 35); they are approximate for other settings.  If adjacent epochs overlap after the correction, the overlap is split at its midpoint.

\bold{Tuning and accuracy.}  The defaults for \code{method="shift"} were tuned on dyads simulated by \code{SymmetryEpochTest-261111.R}: 1000-sample series in which epochs of 10 to 200 samples, separated by gaps of 0 to 50 samples of independent noise, contain lagged copies of a noisy sine, sines of different periods, or a weak lagged regression.  With wMax = 30 and tMax = 30, about three quarters of the epochs longer than 60 samples were found, with start and end typically within 10 samples.  Epochs shorter than about wMax + tMax samples are seldom found, because no window lies wholly inside them.

A fixed threshold has no universal value because the dispersion depends on the noise in the data, the autocorrelation of the series, the window size, and the number of lags.  The default of .6 suits noisy data in which non-symmetric intervals contain noise in both series.  When one series stays a smooth oscillation while the other is disrupted, the cross correlations across lags stay smooth and oscillating, a shift in lag can still align them reasonably well, and the dispersion of disrupted windows may stay below .6.  For such clean signals use \code{threshold="auto"}.  A threshold with a probabilistic interpretation can be taken from surrogate data: pair one series with a time-shifted or unrelated version of the other, calculate the dispersion for the surrogate pairs in the same way, and use a low quantile (for example the 5th percentile) of the surrogate dispersions as the threshold.  The examples show this.

Symmetry in this sense does not require high correlation.  A weak cross correlation pattern that holds still is also symmetric.  Inspect the correlations within an epoch, or test them against surrogates with \code{\link{wccTestSurrogateDyads}}, before interpreting an epoch as an interval of association.
}

\value{
A list with two elements:

\code{epochs}, a data frame with one row per epoch of symmetry and the following columns:
\tabular{ll}{
\sQuote{formationWindow}  \tab  The window at which symmetry forms (the first symmetric window).\cr

\sQuote{breakingWindow}  \tab  The window at which symmetry breaks (the first window after the epoch).  NA if the epoch lasts to the last window.\cr

\sQuote{formationTime}  \tab  The elapsed time of \code{formationWindow}.\cr

\sQuote{breakingTime}  \tab  The elapsed time of \code{breakingWindow}.  NA if the epoch lasts to the last window.\cr

\sQuote{duration}  \tab  The number of windows in the epoch.\cr

\sQuote{startTime}  \tab  The estimated time, in the time series, at which symmetry forms.  Because a window's elapsed time is the end of its data, symmetry is first seen in the windows some time after it forms in the series; \code{startTime} subtracts that delay from \code{formationTime}.  See Details.\cr

\sQuote{endTime}  \tab  The estimated time, in the time series, at which symmetry breaks, or the last elapsed time if the epoch lasts to the end of the series.\cr

\sQuote{meanDispersion}  \tab  The mean (unsmoothed) dispersion over the windows of the epoch.\cr

\sQuote{medianDriftRate}  \tab  The median drift rate over the windows of the epoch.  Near 0 for a stationary epoch.\cr

\sQuote{peakCorrelation}  \tab  The median over the windows of the epoch of the largest cross correlation in each window.\cr

\sQuote{peakLag}  \tab  The median lag of that peak.\cr

\sQuote{peakLagStability}  \tab  The proportion of windows in the epoch whose peak is within one lag column of \code{peakLag}.  Near 1 for a peak that holds still; lower for a pattern that drifts or moves.\cr

\sQuote{peakSharpness}  \tab  The median sharpness of the peak (see \code{\link{wccZeroChangeDirections}}).  A sharp peak indicates a component shared by the two series at that lag that is not autocorrelated.\cr
}
The data frame has zero rows if no epoch is found.

\code{windows}, the data frame returned by \code{wccZeroChangeDirections} with two added columns: \code{smoothedDispersion}, the dispersion after the running median, and \code{symmetric}, a logical indicating whether each window belongs to an epoch.

\code{threshold}, the threshold that was used.
}

\references{
Boker, S. M., Rotondo, J. L., Xu, M., & King, K. (2002). Windowed cross-correlation and peak picking for the analysis of variability in the association between behavioral time series. Psychological methods, 7(3), 338.}

\seealso{
\code{\link{wccZeroChangeDirections}} for the calculation of the dispersion of zero-change directions.  \code{\link{wccCalc}} for the definition of the windowed cross correlation matrix.  \code{\link{wccTestSurrogateDyads}} for surrogate tests of association.
}

\examples{
# Noisy dyad with three epochs of symmetry separated by independent noise:
# a lagged copy of a noisy sine (samples 101-250), sines of different
# periods whose relative phase drifts (301-480), and a weak lagged
# regression (531-700).  The weak regression is the hardest to detect and
# may be reported in pieces.
set.seed(3)
n <- 1000
x <- rnorm(n)
y <- rnorm(n)
s <- sin(1:n/10) + rnorm(n, sd=.5)
x[101:250] <- s[101:250];  y[101:250] <- s[(101:250) + 8]
s1 <- sin(1:n/9) + rnorm(n, sd=.5);  s2 <- sin(1:n/11) + rnorm(n, sd=.5)
x[301:480] <- s1[301:480]; y[301:480] <- s2[(301:480) + 5]
s3 <- sin(1:n/10) + rnorm(n, sd=.5); r3 <- .5 * s3 + rnorm(n, sd=.5)
x[531:700] <- s3[531:700]; y[531:700] <- r3[(531:700) - 6]

wccOut <- wccCalc(x, y, wMax=30, tMax=30)
sym <- wccSymmetryEpochs(wccOut, wMax=30, tMax=30)
sym$epochs[, c("startTime", "endTime", "medianDriftRate", "peakCorrelation",
               "peakLagStability", "peakSharpness")]

# Clean signals: two sines of slightly different frequency, whose phase
# relation drifts uniformly (a symmetry), with two intervals in which
# the second series is replaced by noise (symmetry breaking).  For clean
# signals the threshold is better chosen from the data.
set.seed(1)
tSeries1 <- sin(c(1:1000)/10) + rnorm(1000, mean=0, sd=.005)
tSeries2 <- sin(c(11:1010)/11) + rnorm(1000, mean=0, sd=.005)
tSeries2[201:300] <- rnorm(100, mean=0, sd=1)
tSeries2[501:600] <- rnorm(100, mean=0, sd=1)
wccOut <- wccCalc(tSeries1, tSeries2, wMax=30, tMax=50, wInc=1, tInc=1)
sym <- wccSymmetryEpochs(wccOut, wMax=30, tMax=50, threshold="auto")
sym$threshold
sym$epochs[, c("startTime", "endTime", "medianDriftRate")]

# Split epochs wherever the drift rate changes by more than .02 lag per sample.
wccSymmetryEpochs(wccOut, wMax=30, tMax=50, threshold="auto",
     driftTolerance=.02)$epochs[, c("startTime", "endTime", "medianDriftRate")]

# A threshold from surrogates: pair x with circularly shifted copies of y,
# which keeps each series' autocorrelation but removes their coupling.
surrogateDispersion <- unlist(lapply(1:5, function(k) {
    yShift <- c(y[-(1:(150*k))], y[1:(150*k)])
    wccSurr <- wccCalc(x, yShift, wMax=30, tMax=30)
    wccZeroChangeDirections(wccSurr, wMax=30, tMax=30)$dispersion
}))
surrogateThreshold <- quantile(surrogateDispersion, .10, na.rm=TRUE)
surrogateThreshold
wccSymmetryEpochs(wccCalc(x, y, wMax=30, tMax=30), wMax=30, tMax=30,
     threshold=surrogateThreshold)$epochs[, c("startTime", "endTime")]
}

\keyword{ts}
