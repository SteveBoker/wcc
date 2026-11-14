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

\name{wccCalc}
\alias{wccCalc}

\title{Windowed Cross Correlation function}

\description{
This function calculates a windowed cross correlation matrix from two time series
}

\usage{
wccCalc(inSeries1, inSeries2, wMax=50, tMax=50, wInc=1, tInc=1,
        method=c("c", "cumr", "cumc", "cumcuda", "r"),
        precision=c("double", "single"), ...)
}

\arguments{
\item{inSeries1}{A vector of numeric values. This is the first time series on which windowed cross correlation is calculated.  When inSeries1 is leading \code{inSeries2}, the peak correlation closest to zero in negative in wccPlot.  Note that \code{inSeries1} must be of the same length as \code{inSeries2}.   Note that the only difference between \code{inSeries1} and \code{inSeries2} is which is treated as positive lag and which is negative lag in wccPlot.}
\item{inSeries2}{A vector of numeric values. This is the second time series on which windowed cross correlation is calculated. Must be of the same length as \code{inSeries1}.   When \code{inSeries2} is leading \code{inSeries1}, the peak correlation closest to zero in positive in wccPlot. }
\item{wMax}{A numeric indicating the number of samples in a window (default=50).  This must be greater or equal to 5 samples.  It is recommended that wMax be greater than 15.}
\item{tMax}{The maximum number of samples to lag the windows (default=50).  The cross correlation is calculated for windows that are time lagged against each other, allowing for processes than require some time to evolve.  For instance, if the two time series are from motion capture of two individuals conversing with one another, the speaker's movements are likely to occur prior to the responses of the listener.  tMax should be set to be greater than the number of samples that represents the maximum time lag that is likely to occur. }
\item{wInc}{The number of samples incremented between successive windows  (default=1).  This is most often set to 1.  If the sample rate is very high and the time series is long, wInc may be increased.  However, \code{wccCalc}'s maximum sensitivity to nonstationary change in the time series occurs when wInc is 1.}
\item{tInc}{The number of samples incremented between successive lags  (default=1).  This is most often set to 1.  If the sample rate is high and the maximum lag is high then tInc may be increased.  However, \code{wccCalc}'s maximum sensitivity to changes in lags occurs when the tInc is set to 1.}
\item{method}{Backend selector for the windowed cross correlation.  \code{"c"} (default) calls the original compiled \code{windcross} binary; \code{"r"} is the original pure-R loop.  \code{"cumc"} uses a cumulative-sum C backend that computes each cell in O(1) from prefix sums, giving O(n * nLags) cost independent of window size; \code{"cumr"} is the pure-R reference for that algorithm.  \code{"cumcuda"} runs the same cumulative-sum algorithm on the GPU and requires an installation built with CUDA support (see \code{src/Makevars}); without it the call raises an informative error.  The \code{cum*} methods do not support missing data: any NA in the input is an error.  When \code{tInc > 1} the \code{cum*} methods use lags \code{0, tInc, 2*tInc, ..., tMax}, matching \code{"c"}; the legacy \code{"r"} method shifts by the column index instead.}
\item{precision}{Floating-point precision of the GPU computation; only valid with \code{method="cumcuda"}.  \code{"double"} (default) computes in FP64.  \code{"single"} computes in FP32, which is substantially faster on consumer GPUs (whose FP64 throughput is a small fraction of FP32) at the cost of correlations accurate to roughly 1e-5 rather than machine precision.  Results are always returned as R doubles.}
\item{...}{Reserved for the deprecated \code{windcross=TRUE/FALSE} argument, mapped to \code{method="c"} / \code{method="r"} with a warning.}
}

\details{
\code{wccCalc} calculates cross correlations between windows, i.e., sampled sections of two time series.  These windows are sampled at time lagged intervals from one another so that the association between the time series can be estimated even though the time lag of maximum association between the two time series may be itself time varying.  For instance, when two individuals converse, their body motions may be associated such that the leader and follower take turns as speaker and listener.  This same kind of nonstationary association has been observed between many physiological time series.  At any particular elapsed time, the maximum correlation closest to a lag of zero (whether a positive or negative lag) is taken to be the lagged offset between the two time series.  For series that have oscillatory structure (such as head nodding in conversation)  the maximum correlation closest to a lag of zero can be interpreted as the phase lag between the series.

Note that missing values in \code{inSeries1} and \code{inSeries2} are treated such that each window calculates correlation on pairwise complete observations.  If there are fewer than 5 non-missing pairs in any window then the correlation for that window is returned as NA.
}
    
\value{
Returns an N x P matrix of numeric values of cross correlations, where N is the number of windows calculated and P is 2*floor(tMax/tInc) + 1, the total lagged values of each window.  Each row of the returned matrix represents one elapsed time at sample T.  The number of columns is always odd and the center column, C = floor(tMax/tInc) + 1 is the cross correlation for zero lag for the windows from \code{inSeries1} and \code{inSeries2} whose last element occurs at elapsed sample T.  Columns L less than the center column contain the cross correlation values for the windows where the window from \code{inSeries1} occurs \code{L*tInc} samples earlier than the window from \code{inSeries2} whose last value is at sample T.  Similarly columns L greater than the center column contain cross corrlations for windows where the window from \code{inSeries2} occurs \code{L*tInc} samples earlier than the window from \code{inSeries1} whose last value is at sample T.
}

\references{
Boker, S. M., Rotondo, J. L., Xu, M., & King, K. (2002). Windowed cross-correlation and peak picking for the analysis of variability in the association between behavioral time series. Psychological methods, 7(3), 338.}


\examples{
#Calculate windowed cross correlation between two time series
tSeries1 <- sin(c(1:1000)/10) + rnorm(1000, mean=0, sd=.5)
tSeries2 <- sin(c(1:1000)/15) + rnorm(1000, mean=0, sd=.5)
wccCalcOut <- wccCalc(inSeries1=tSeries1, inSeries2=tSeries2, wMax=50, tMax=50, wInc=1,
     tInc=1, method="r")

}
