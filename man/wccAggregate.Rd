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

\name{wccAggregate}
\alias{wccAggregate}

\title{wccAggregate}

\description{
This function runs \code{wccCalc} and peakpicking on a pair of time series and returns aggregate statistics.
}

\usage{
wccAggregate(inSeries1=NA, inSeries2=NA, wMax=50, tMax=50, wInc=1, tInc=1, Lsize=8,
     pspan=.25, type="Max", samplespersecond=1,
     method=c("c", "cumr", "cumc", "cumcuda", "r"), embedD=9, ...)
}

\arguments{
\item{inSeries1}{A vector of numeric values. This is the first time series on which windowed cross correlation is calculated.  When \code{inSeries1} is leading \code{inSeries2}, the peak correlation closest to zero in negative in \code{wccPlot}.  Note that \code{inSeries1} must be of the same length as \code{inSeries2}.   Note that the only difference between \code{inSeries1} and \code{inSeries2} is which is treated as positive lag and which is negative lag in \code{wccPlot}.}
\item{inSeries2}{A vector of numeric values. This is the second time series on which windowed cross correlation is calculated.  Must be of the same length as \code{inSeries1}.   When \code{inSeries2} is leading \code{inSeries1}, the peak correlation closest to zero in positive in \code{wccPlot}.}
\item{wMax}{A numeric indicating the number of samples in a window (default=50).  This must be greater or equal to 5 samples.  It is recommended that wMax be greater than 15.}
\item{tMax}{The maximum number of samples to lag the windows (default=50).  The cross correlation is calculated for windows that are time lagged against each other, allowing for processes than require some time to evolve.  For instance, if the two time series are from motion capture of two individuals conversing with one another, the speaker's movements are likely to occur prior to the responses of the listener.  tMax should be set to be greater than the number of samples that represents the maximum time lag that is likely to occur. }
\item{wInc}{The number of samples incremented between successive windows (default=1).  This is most often set to 1.  If the sample rate is very high and the time series is long, wInc may be increased.  However, \code{wccCalc}'s maximum sensitivity to nonstationary change in the time series occurs when wInc is 1.}
\item{tInc}{The number of samples incremented between successive lags  (default=1).  This is most often set to 1.  If the sample rate is high and the maximum lag is high then tInc may be increased.  However, \code{wccCalc}'s maximum sensitivity to changes in lags occurs when the tInc is set to 1.}
\item{Lsize}{An integer specifying the width required for a peak to be identified by \code{wccPeakPick} (default=8).  An Lsize that is too small may find spurious peaks that are just noise, whereas an Lsize that is too large may miss actual peaks.}
\item{pspan}{A numeric specifying the amount of smoothing applied by the \code{\link{loess}} function within \code{wccPeakPick}  (default=.25).}
\item{type}{A character string indicating whether a maximum, \sQuote{Max}, or minimum, \sQuote{Min}, should be detected by \code{wccPeakPick} (default="Max").}
\item{samplespersecond}{A numeric indicating the sampling rate of the time series.  This will determine the units of the first derivative aggregate statistics. (default=1)}
\item{method}{Backend passed through to \code{\link{wccCalc}}: \code{"c"} (default, original compiled windcross), \code{"r"} (original pure-R loop), \code{"cumc"} (cumulative-sum C backend), \code{"cumr"} (cumulative-sum pure-R reference), or \code{"cumcuda"} (cumulative-sum CUDA backend, requires a CUDA-enabled build).  See \code{\link{wccCalc}} for details and the missing-data caveat for the \code{cum*} methods.}
\item{embedD}{A numeric indicating the embedding dimension passed to GLLAEmbed to calculate derivatives of the peak picking lag timeseries.  The default (9) will provide a medium smooth at the cost of some loss of resolution in the timing of jumps in the second derivative.  The smallest legal value is 5, which provides less smoothing and the best time resolution.  If your sampling rate is low, then a smaller value of embedD is useful, while if your sampling rate is high, a larger value of embedD will provide better noise rejection by smoothing.}
\item{...}{Reserved for the deprecated \code{windcross=TRUE/FALSE} argument, mapped to \code{method="c"} / \code{method="r"} with a warning.}
}

\details{
wccAggregate is a wrapper for the \code{wccCalc} and wccPeakpick functions.  The returns from these two function calls are then aggregated into 10 different statistics representing characteristics of \code{wccCalc} and wccPeakpick for the pair of input time series.  wccAggregate also must be called with a filename in the \code{savefile} argument prior to calling \code{wccPlot} so that the raw results can be plotted.  The wccAggregate function returns a single row dataframe so that it can be called repeatedly by the wccSurrogateDyad function to create a distribution of surrogate Windowed Cross Correlation results.  The arguments are also saved in this dataframe so that the wccFindDyadParam function can explore the space of arguments to find appropriate values from a hold-out set of timeseries to then be used to calculate surrogate distributions for the full data set.
}
    
\value{
wccAggregate returns a data frame with one row and the following columns:
\tabular{ll}{
\sQuote{wMax}   \tab  The value of wMax in the argument passed to \code{wccCalc}.\cr

\sQuote{tMax}   \tab  The value of tMax in the argument passed to \code{wccCalc}. \cr

\sQuote{wInc}   \tab  The value of wInc in the argument passed to \code{wccCalc}. \cr

\sQuote{tInc}   \tab   The value of tInc in the argument passed to \code{wccCalc}.\cr

\sQuote{Lsize}   \tab  The value of Lsize in the argument passed to \code{wccPeakPick}. \cr

\sQuote{pspan}   \tab  The value of pspan in the argument passed to \code{wccPeakPick}. \cr

\sQuote{type}   \tab   The value of type as a numeric translation of the type argument passed to \code{wccPeakPick}.  This is +1 when type = \sQuote{Max} and -1 when type = \sQuote{Min}.\cr

\sQuote{samples}   \tab  Length of the time series \code{inSeries1} and \code{inSeries2}\cr

\sQuote{windows}   \tab  The number of windows returned by \code{wccCalc}.\cr

\sQuote{pctMissing}   \tab  Proportion of missing values in the matrix returned by \code{wccCalc}.\cr

\sQuote{pctMissingWindows}   \tab  Proportion of windows (i.e., rows) in the matrix returned by \code{wccCalc} that include at least one missing value.\cr

\sQuote{maxMean}   \tab  Mean of the Fischer's Z transform of the maximum correlation vector returned returned by \code{wccPeakPick}.  Note that if type=\sQuote{Min} then this value is the Fischer's Z transform of the minimum correlation vector. \cr

\sQuote{maxVar}   \tab   Variance of the Fischer's Z transform of the maximum correlation vector returned returned by \code{wccPeakPick}.  Note that if type=\sQuote{Min} then this value is the Fischer's Z transform of the minimum correlation vector. \cr

\sQuote{totalMean}   \tab  Mean of the Fischer's Z transform of all values in the matrix returned by \code{wccCalc}.\cr

\sQuote{totalVar}   \tab  Variance of the Fischer's Z transform of all values in the matrix returned by \code{wccCalc} .\cr

\sQuote{zeroLagMean}   \tab  Mean of the Fischer's Z transform of all values in the zero lag column of the matrix returned by \code{wccCalc}.\cr

\sQuote{zeroLagVar}   \tab Variance of the Fischer's Z transform of all values in the zero lag column of the matrix returned by \code{wccCalc}. \cr

\sQuote{lagMean}   \tab  Mean of the lag vector returned by \code{wccPeakPick}.\cr

\sQuote{lagVar}   \tab  Variance of the lag vector returned by \code{wccPeakPick}.\cr

\sQuote{dlagMean}   \tab Mean of the first derivative of the lag vector returned by \code{wccPeakPick}. \cr

\sQuote{dlagVar}   \tab  Variance of the first derivative of the lag vector returned by \code{wccPeakPick}.\cr

\sQuote{d2lagMean}   \tab Mean of the second derivative of the lag vector returned by \code{wccPeakPick}. \cr

\sQuote{d2lagVar}   \tab  Variance of the second derivative of the lag vector returned by \code{wccPeakPick}.\cr

\sQuote{elapsedSeconds}   \tab  Total elapsed seconds required to perform all calculations and function calls.\cr

}
}

\references{
Boker, S. M., Rotondo, J. L., Xu, M., & King, K. (2002). Windowed cross-correlation and peak picking for the analysis of variability in the association between behavioral time series. Psychological methods, 7(3), 338.}

\seealso{
\code{\link{wccCalc}}, \code{\link{wccPeakPick}}, \code{\link{loess}}. 
}

\examples{
#Calculate aggregated statistics for windowed cross correlation between two time series
tSeries1 <- sin(c(1:1000)/10) + rnorm(1000, mean=0, sd=.5)
tSeries2 <- sin(c(1:1000)/15) + rnorm(1000, mean=0, sd=.5)
wccAggregate(inSeries1=tSeries1, inSeries2=tSeries2, wMax=50, tMax=50, wInc=1, tInc=1,
     Lsize=8, pspan=.25, type="Max", samplespersecond=10, method="c")

}
