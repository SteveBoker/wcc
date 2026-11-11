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

\name{wccZeroChangeDirections}
\alias{wccZeroChangeDirections}

\title{wccZeroChangeDirections}

\description{
This function estimates, for each window of a windowed cross correlation matrix, the direction in the elapsed time by lag plane along which the cross correlation does not change (the isobar), and how well a single such direction describes the change at every lag.  Low disagreement indicates symmetry: the pattern of cross correlations across lags is either stationary or drifting in lag at a uniform rate.  The result is used by \code{\link{wccSymmetryEpochs}} to estimate when symmetry forms and when it breaks.
}

\usage{
wccZeroChangeDirections(wccMatrix, wMax=NA, tMax=NA, wInc=1, tInc=1,
     samplespersecond=1, method=c("shift", "gradient"), timeScale=NA,
     maxShift=NA, timeSmooth=NA, normalizeRows=TRUE,
     weight=c("gradient", "none"), minGrad=1e-8)
}

\arguments{
\item{wccMatrix}{A numeric matrix of windowed cross correlations as returned by \code{\link{wccCalc}}: rows are windows (elapsed time) and columns are lags.  It must have at least 3 rows and 3 columns.}
\item{wMax}{The number of samples in a window that was used in \code{wccCalc} (default=NA).  This is used to calculate the elapsed time of each window and the default values of \code{timeScale} and \code{timeSmooth}.  It is required for \code{method="shift"} unless \code{timeScale} is given.  If NA, elapsed time starts at 0.}
\item{tMax}{The maximum lag that was used in \code{wccCalc} (default=NA).  This is used to calculate the elapsed time of each window.}
\item{wInc}{The number of samples between successive windows that was used in \code{wccCalc} (default=1).}
\item{tInc}{The number of samples between successive lags that was used in \code{wccCalc} (default=1).}
\item{samplespersecond}{A numeric indicating the sampling rate of the time series (default=1).  Elapsed times are returned in seconds.  With the default of 1, elapsed times are in samples.}
\item{method}{Either \sQuote{shift} (default) or \sQuote{gradient}.  See Details.  \sQuote{shift} compares the cross correlations across lags one window length apart and finds the single shift in lag that best explains the change.  \sQuote{gradient} uses the local gradient at every cell.  \sQuote{shift} detected the start and end of known epochs of symmetry in simulated data much more accurately than \sQuote{gradient}.}
\item{timeScale}{For \code{method="shift"}: the number of windows, h, before and after each window that are compared (default=NA).  When NA, \code{wMax/(2*wInc)} is used, so that the two windows compared (h windows before and h windows after) do not share data.  Ignored for \code{method="gradient"}.}
\item{maxShift}{For \code{method="shift"}: the largest shift in lag, in lag columns, that is tried when the window h windows ahead is compared with the window h windows behind (default=NA).  When NA, \code{ceiling(timeScale/2)} is used, which allows drift of up to a quarter of a lag column per window.  It bounds the drift rate that can be recognized as symmetry, to about \code{maxShift / (2 * timeScale)} lag columns per window.  Ignored for \code{method="gradient"}.}
\item{timeSmooth}{The width, in windows, of a centered moving average applied along elapsed time before the directions are calculated (default=NA).  For \code{method="shift"}, NA means 1 (no smoothing).  For \code{method="gradient"}, NA means \code{wMax/(4*wInc)} (a quarter of a window length) is used if \code{wMax} is given, otherwise 1.  Successive windows overlap and each is a noisy estimate of the cross correlation, so without smoothing the window-to-window sampling noise is read as change.  How much smoothing is needed depends on how noisy the cross correlations are.  For clean, nearly deterministic signals little or no smoothing gives the sharpest boundaries.  For noisy behavioral data, smoothing over up to one window length (\code{wMax/wInc}) separates symmetric from non-symmetric windows better, at the cost of blurring the formation and breaking times.  Set to 1 for no smoothing.}
\item{normalizeRows}{A logical (default=TRUE).  If TRUE, each window is standardized across lags (mean 0, standard deviation 1) before the gradients are calculated.  This makes the measure depend only on the shape and position of the cross correlation pattern across lags, not on its height.  Without this, a peak that stays at the same lag but rises and falls in height is counted as broken symmetry.  Standardization is done before the smoothing given by \code{timeSmooth}, so that windows with low cross correlations, for example where one series is pure noise, are not swamped by neighboring windows with high cross correlations.  Setting this to FALSE is rarely useful with real data.}
\item{weight}{A character string, either \sQuote{gradient} (default) or \sQuote{none}.  With \sQuote{gradient}, each cell is weighted by its squared gradient magnitude, so that cells where the surface is flat, such as the crest of a peak, where the direction of zero change is undefined, contribute little.  With \sQuote{none}, all cells are weighted equally and cells whose gradient magnitude is below \code{minGrad} are dropped.}
\item{minGrad}{A numeric (default=1e-8).  Cells with gradient magnitude below this value are dropped when \code{weight="none"}.  Ignored when \code{weight="gradient"}.}
}

\details{
If the cross correlation pattern across lags is translating rigidly in lag at a rate v, so that r(t, lag) = f(lag - v t), then the cross correlation does not change along lines of slope dlag/dt = v at every lag.  A single direction of zero change describes the whole pattern, and the pattern is symmetric.  A stationary pattern is the special case v = 0.  This definition of symmetry therefore includes peaks that drift in lag at a uniform rate, as when one person's lead over another grows steadily.  When the pattern changes in any other way, for example when one peak fades while another grows, no single direction describes it.

\bold{method="shift"} (default).  Each window's cross correlations are first standardized across lags.  For window t, the pattern h = \code{timeScale} windows later is compared with the pattern h windows earlier.  The earlier pattern is shifted in lag by every whole number of lag columns from \code{-maxShift} to \code{maxShift}, and the shift with the highest correlation between the two patterns is kept.  That shift, divided by 2h, is the slope of the direction of zero change.  The \code{coherence} is the correlation at that shift: 1 when the later pattern is exactly the earlier one translated uniformly in lag.  The \code{dispersion} is 1 - coherence.  Because windows of \code{wccCalc} overlap, neighboring windows share most of their data and are nearly identical even when the two time series are unrelated noise.  With the default h = wMax/(2 wInc), the two windows compared share no data, so noise is not mistaken for a stationary pattern.  The first and last h windows have no score (NA).

\bold{method="gradient"}.  For each cell of the matrix the gradient of the cross correlation, (dr/dt, dr/dlag), is calculated by central differences.  The direction perpendicular to the gradient is the direction in which the cross correlation does not change.  Its slope, dlag/dt = -(dr/dt)/(dr/dlag), is the rate at which a contour of the surface moves in lag as elapsed time increases.

For a rigidly translating pattern this slope equals v at every lag, so all the directions within the window agree and their dispersion is 0.  When the pattern changes shape the directions disagree and the dispersion approaches 1.  The gradient method compares neighboring windows, which share most of their data; it works well for clean signals but can mistake unrelated noise for a stationary pattern.

For the gradient method, directions of zero change are axial (a line has no arrow), so the angles are doubled before they are averaged.  With the default gradient weighting, the dispersion equals one minus the coherence of the window's 2 by 2 structure tensor.

Both methods work in units of windows and lag columns, not in seconds.  This keeps the dispersion independent of the choice of units for elapsed time and lag.  The arguments \code{wInc}, \code{tInc}, and \code{samplespersecond} only rescale the reported elapsed times and drift rates.

The direction of the drift follows the column order of \code{wccCalc}.  A positive \code{driftRate} means the pattern is moving toward lags at which \code{inSeries2} leads \code{inSeries1} by more, or lags by less.  A negative \code{driftRate} means it is moving toward lags at which \code{inSeries1} leads \code{inSeries2} by more.

Symmetry in this sense does not require high correlation.  A weak cross correlation pattern that holds still is also symmetric.  Use \code{peakCorrelation}, or surrogate tests to decide whether a symmetric interval is also an interval of meaningful association.
}

\value{
A data frame with one row per window of \code{wccMatrix} and the following columns:
\tabular{ll}{
\sQuote{window}  \tab  The row index of the window in \code{wccMatrix}.\cr

\sQuote{elapsedTime}  \tab  The elapsed time of the end of the window, (wMax + tMax + (window - 1) * wInc) / samplespersecond.  If \code{wMax} or \code{tMax} is NA, the offset wMax + tMax is taken as 0.\cr

\sQuote{dispersion}  \tab  How poorly a single direction of zero change describes the change at every lag: 0 for symmetry.  For \code{method="shift"}, 1 minus the best correlation between the patterns h windows apart, from 0 to 2; NA for the first and last h windows.  For \code{method="gradient"}, the dispersion of the local directions, from 0 to 1; NA when the window has no usable gradient.\cr

\sQuote{coherence}  \tab  1 - dispersion.\cr

\sQuote{direction}  \tab  The direction of zero change in radians, measured from the elapsed time axis toward increasing lag column, between -pi/2 and pi/2.\cr

\sQuote{driftRate}  \tab  The drift in lag per unit of elapsed time implied by \code{direction}, tan(direction) * tInc / wInc.  Values near 0 indicate a stationary pattern.  The value is only meaningful when \code{dispersion} is low.\cr

\sQuote{gradEnergy}  \tab  For \code{method="gradient"}, the mean squared gradient magnitude in the window, after any smoothing and normalization.  NA for \code{method="shift"}.\cr

\sQuote{peakCorrelation}  \tab  The largest cross correlation in the window.\cr

\sQuote{peakLag}  \tab  The lag of \code{peakCorrelation}, in lag units (tInc / samplespersecond per column), negative toward the first columns of \code{wccCalc}.\cr

\sQuote{peakSharpness}  \tab  \code{peakCorrelation} minus the mean of the cross correlations 3 lag columns either side of the peak.  A sharp peak indicates a component shared by the two series at that lag that is not autocorrelated.\cr
}
The parameters used are stored in the attribute \code{"wccSymmetryParams"}, which \code{\link{wccSymmetryEpochs}} uses to convert formation and breaking times to estimated start and end times.
}

\references{
Boker, S. M., Rotondo, J. L., Xu, M., & King, K. (2002). Windowed cross-correlation and peak picking for the analysis of variability in the association between behavioral time series. Psychological methods, 7(3), 338.}

\seealso{
\code{\link{wccSymmetryEpochs}} to find the elapsed times of symmetry formation and symmetry breaking.  \code{\link{wccCalc}} for the definition of the windowed cross correlation matrix.  \code{\link{wccPeakPick}} for tracking the lag of a single peak.
}

\examples{
# Two time series that are uncoupled, then coupled at a fixed lag,
# then coupled at a lag that grows steadily.
set.seed(2)
N <- 3000
x <- as.numeric(arima.sim(list(ar=.9), N))
y <- as.numeric(arima.sim(list(ar=.9), N))
for (i in 800:1500) y[i] <- x[i - 6] + .3 * rnorm(1)
for (i in 1800:2400) y[i] <- x[i - round((i - 1800) / 40)] + .3 * rnorm(1)

wccOut <- wccCalc(x, y, wMax=100, tMax=20, wInc=10, tInc=1)
zcd <- wccZeroChangeDirections(wccOut, wMax=100, tMax=20, wInc=10)

# Dispersion is low while the pattern is stationary or drifting uniformly.
plot(zcd$elapsedTime, zcd$dispersion, type='l',
     xlab="Elapsed Time (samples)", ylab="Dispersion of zero-change directions")
abline(v=c(800, 1500, 1800, 2400), lty=2)

# The gradient method on the same matrix, for comparison.
zcdG <- wccZeroChangeDirections(wccOut, wMax=100, tMax=20, wInc=10,
     method="gradient", timeSmooth=10)
lines(zcdG$elapsedTime, zcdG$dispersion, col="red")
}

\keyword{ts}
