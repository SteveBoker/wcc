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

\name{wccPlot}
\alias{wccPlot}

\title{wccPlot}

\description{
This function plots a Windowed Cross Correlation heatmap.
}

\usage{
wccPlot(inSeries1=NA, inSeries2=NA, startwindow=1, endwindow=200, wMax=50, tMax=50,
    wInc=1, tInc=1, Lsize=8, pspan=.25, type="Max", samplespersecond=1,
    method=c("c", "cumr", "cumc", "cumcuda", "r"), ...)
}

\arguments{
\item{inSeries1}{A vector of numeric values. This is the first time series on which windowed cross correlation is calculated.  When \code{inSeries1} is leading \code{inSeries2}, the peak correlation closest to zero in negative in \code{wccPlot}.  Note that \code{inSeries1} must be of the same length as \code{inSeries2}.   Note that the only difference between \code{inSeries1} and \code{inSeries2} is which is treated as positive lag and which is negative lag in \code{wccPlot}.}
\item{inSeries2}{A vector of numeric values. This is the second time series on which windowed cross correlation is calculated. Must be of the same length as \code{inSeries1}.   When \code{inSeries2} is leading \code{inSeries1}, the peak correlation closest to zero in positive in \code{wccPlot}. }
\item{startwindow}{A numeric value specifying the first window to plot (default=1).  Windows are plotted as vertical stripes in the heatmap and the elapsed time of this window is \code{wMax + tMax + (startwindow * wInc / samplespersecond)}}
\item{endwindow}{A numeric value specifying the last window to plot (default=200).}
\item{wMax}{A numeric indicating the number of samples in a window (default=50).  This must be greater or equal to 5 samples.  It is recommended that wMax be greater than 15.}
\item{tMax}{The maximum number of samples to lag the windows (default=50).  The cross correlation is calculated for windows that are time lagged against each other, allowing for processes than require some time to evolve.  For instance, if the two time series are from motion capture of two individuals conversing with one another, the speaker's movements are likely to occur prior to the responses of the listener.  tMax should be set to be greater than the number of samples that represents the maximum time lag that is likely to occur. }
\item{wInc}{The number of samples incremented between successive windows  (default=1).  This is most often set to 1.  If the sample rate is very high and the time series is long, wInc may be increased.  However, \code{wccCalc}'s maximum sensitivity to nonstationary change in the time series occurs when wInc is 1.}
\item{tInc}{The number of samples incremented between successive lags  (default=1).  This is most often set to 1.  If the sample rate is high and the maximum lag is high then tInc may be increased.  However, \code{wccCalc}'s maximum sensitivity to changes in lags occurs when the tInc is set to 1.}
\item{Lsize}{An integer specifying the width required for a peak to be identified by \code{wccPeakPick} (default=8).  An Lsize that is too small may find spurious peaks that are just noise, whereas an Lsize that is too large may miss actual peaks.}
\item{pspan}{A numeric specifying the amount of smoothing applied by the \link{loess} function within \code{wccPeakPick}  (default=.25).}
\item{type}{A character string indicating whether a maximum, \sQuote{Max}, or minimum, \sQuote{Min}, should be detected by \code{wccPeakPick} (default=\sQuote{Max}). If type=\sQuote{Zero} then a topographic map is generated rather than a line from wccPeakPick\code{wccPeakPick}.}
\item{samplespersecond}{A numeric indicating the sampling rate of the time series.  This will determine the units of the first derivative aggregate statistics. (default=1)}
\item{method}{Backend passed through to \code{\link{wccCalc}}: \code{"c"} (default, original compiled windcross), \code{"r"} (original pure-R loop), \code{"cumc"} (cumulative-sum C backend), \code{"cumr"} (cumulative-sum pure-R reference), or \code{"cumcuda"} (cumulative-sum CUDA backend, requires a CUDA-enabled build).  See \code{\link{wccCalc}} for details and the missing-data caveat for the \code{cum*} methods.}
\item{...}{Reserved for the deprecated \code{windcross=TRUE/FALSE} argument, mapped to \code{method="c"} / \code{method="r"} with a warning.}
}

\details{
This function runs \code{wccCalc} and \code{wccPeakPick} using the arguments passed in and plots a heatmap of a selected section of windows.  The vertical axis is labeled with the lags in the units of seconds.  Each calculated window is plotted as a vertical slice of the heat map and correlations near 1 are plotted in yellow and correlations near -1 are plotted in red.  Lags of \code{inSeries1} are plotted as positive and lags of \code{inSeries2} are plotted as negative lags.    The horizontal black line in the middle of the plot represents a lag of zero.    Thus, if a feature in \code{inSeries1} occurs prior to a similar feature in \code{inSeries2}, there will be a positive correlation shown at a negative lag, i.e., below the horizontal zero lag line.  The result of \code{wccPeakPick} is overlayed as a black line that traces the peak correlation closest to a lag of zero.  The horizontal axis is labeled with the elapsed time to the right edge of two windows when they are at zero lag.  

It is recommended that the  plot includes no more than 1000 windows since if a large number of windows are plotted, the resulting heatmap will overwrite one window with the next window.  Other considerations include whether the interval of time of interest is near the beginning or end of the timeseries.  The elapsed time associated with the first window will be \code{tMax*tInc*samplespersecond}. Thus an elapsed time of 0 will never occur. It is recommended that the startwindow and endwindow be set to bracket a sample of the timeseries that is of interest to the researcher.  If endwindow - startwindow is about 500, a balance between sensitivity to nonstationarity and plotting resolution is achieved.
}
    
\value{
Returns nothing.
}

\references{
Boker, S. M., Rotondo, J. L., Xu, M., & King, K. (2002). Windowed cross-correlation and peak picking for the analysis of variability in the association between behavioral time series. Psychological methods, 7(3), 338.}

\seealso{
\code{\link{wccCalc}} and \code{\link{wccPeakPick}} for the calculation of the \code{wccCalc} matrix and \code{wccPeakPick} timeseries that are used by this function. 
}

\examples{
#Create a windowed cross correlation plot for two simulated time series with a phase offset
tSeries1 <- sin(c(1:1000)/10) + rnorm(1000, mean=0, sd=.5)
tSeries2 <- sin(1+c(1:1000)/10) + rnorm(1000, mean=0, sd=.5)
wccPlot(inSeries1=tSeries1, inSeries2=tSeries2, startwindow=1, endwindow=500, 
    wMax=100, tMax=100, wInc=1, tInc=1, Lsize=8, pspan=.25, type="Max", 
    samplespersecond=10)
    

}
