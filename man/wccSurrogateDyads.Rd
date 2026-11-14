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

\name{wccSurrogateDyads}
\alias{wccSurrogateDyads}

\title{wccSurrogateDyads}

\description{
This function randomly pairs timeseries and creates a distribution of the null hypothesis that the pairing of dyads does not matter.
}

\usage{
wccSurrogateDyads(inArray1=NA, inArray2=NA, wMax=50, tMax=50, wInc=1, tInc=1, Lsize=8,
     pspan=.25, type="Max", nSurrogates=NA,
     method=c("c", "cumr", "cumc", "cumcuda", "r"), embedD=9, ...)
}

\arguments{
\item{inArray1}{An N by T matrix of numeric values representing N timeseries of length T. These are the first set of time series which are randomly paired with timeseries from \code{inArray2} and on which windowed cross correlation are calculated.  Must be of the same order N by T as \code{inArray2}. Note that the only difference between \code{inArray1} and \code{inArray2} is which is treated as positive lag and which is negative lag in \code{wccPlot}.  When \code{inArray1} is leading \code{inArray2}, the peak correlation closest to zero in negative in \code{wccPlot}.  }
\item{inArray2}{An N by T matrix of numeric values representing N timeseries of length T. These are the second set of time series which are randomly paired with timeseries from \code{inArray1} and on which windowed cross correlation are calculated.  Must be of the same order N by T as \code{inArray1}.   When \code{inArray2} is leading \code{inArray1}, the peak correlation closest to zero in positive in \code{wccPlot}.}
\item{wMax}{A numeric indicating the number of samples in a window (default=50).  This must be greater or equal to 5 samples.  It is recommended that wMax be greater than 15.}
\item{tMax}{The maximum number of samples to lag the windows (default=50).  The cross correlation is calculated for windows that are time lagged against each other, allowing for processes than require some time to evolve.  For instance, if the two time series are from motion capture of two individuals conversing with one another, the speaker's movements are likely to occur prior to the responses of the listener.  tMax should be set to be greater than the number of samples that represents the maximum time lag that is likely to occur. }
\item{wInc}{The number of samples incremented between successive windows (default=1).  This is most often set to 1.  If the sample rate is very high and the time series is long, wInc may be increased.  However, wccCalc's maximum sensitivity to nonstationary change in the time series occurs when wInc is 1.}
\item{tInc}{The number of samples incremented between successive lags  (default=1).  This is most often set to 1.  If the sample rate is high and the maximum lag is high then tInc may be increased.  However, wccCalc's maximum sensitivity to changes in lags occurs when the tInc is set to 1.}
\item{Lsize}{An integer specifying the width required for a peak to be identified by \code{wccPeakPick} (default=8).  An Lsize that is too small may find spurious peaks that are just noise, whereas an Lsize that is too large may miss actual peaks.}
\item{pspan}{A numeric specifying the amount of smoothing applied by the \link{loess} function within \code{wccPeakPick}  (default=.25).}
\item{type}{A character string indicating whether a maximum, \sQuote{Max}, or minimum, \sQuote{Min}, should be detected by \code{wccPeakPick} (default=\sQuote{Max}).}
\item{nSurrogates}{A numeric indicating the number of surrogates to generate (default=NA).  If there are N timeseries in \code{inArray1}, then the maximum number of unique surrogates that can be generated is N*(N-1)/2}
\item{method}{Backend passed through to \code{\link{wccCalc}}: \code{"c"} (default, original compiled windcross), \code{"r"} (original pure-R loop), \code{"cumc"} (cumulative-sum C backend), \code{"cumr"} (cumulative-sum pure-R reference), or \code{"cumcuda"} (cumulative-sum CUDA backend, requires a CUDA-enabled build).  See \code{\link{wccCalc}} for details and the missing-data caveat for the \code{cum*} methods.}
\item{embedD}{A numeric indicating the embedding dimension passed to \code{wccGLLAEmbed} to calculate derivatives of the peak picking lag timeseries.  The default (9) will provide a medium smooth at the cost of some loss of resolution in the timing of jumps in the second derivative.  The smallest legal value is 5, which provides less smooth and the best time resolution.  If your sampling rate is low, then a smaller value of \code{embedD} is useful, while if your sampling rate is high, a larger value of \code{embedD} will provide better noise rejection by smoothing.}
\item{...}{Reserved for the deprecated \code{windcross=TRUE/FALSE} argument, mapped to \code{method="c"} / \code{method="r"} with a warning.}
}

\details{
This function generates and returns a distribution of aggregated statistics from running Windowed Cross Correlation and peakpicking (Boker, Rotondo, Xu, & King 2002) on pairs of timeseries that conform to the null hypothesis that the true pairing of the timeseries does not matter.  This is called a surrogate analysis (Moulder, Boker, Ramsayer, & Tschacher 2018) and is a nonparametric way of estimating whether the observed association between timeseries is due to dyadic interactions occuring while the timeseries were measured or simply due to the distributions of scores within all the timeseries.  This is important to consider since the timeseries may be measurements of individual processes that are constrained in some way.  For instance, consider dyadic timeseries of head motions.  All persons' heads are constrained in the motions they can perform due to similarities in skeletal and muscular dynamics that are shared by most humans.  One does not want these shared constraints to show up as significant associations within dyads when they are actually due to similarities between individuals.

The function returns a dataframe whose columns contain different ways of aggregating the associations and each row contains results from one random pairing of dyads.  The randomization allows reuse of dyads, so it is best if the number of possible random pairings is large relative to the number of surrogates requested.  If there are N rows in  inArray1 and inArray2, the maximum number of unique surrogates that can be generated is N*(N-1)/2.  With that constraint, we recommend at least N = 10 dyads as a source pool since that allows 45 unique surrogates to be generated.   The only exclusion to the random pairing of dyads is that if the random pairing results in an actual dyad, it is not included.
}
    
\value{
Returns a dataframe whose columns are defined in \code{\link{wccAggregate}} and whose rows are the results of each random pairing of dyads.
}

\references{
Boker, S. M., Rotondo, J. L., Xu, M., & King, K. (2002). Windowed cross-correlation and peak picking for the analysis of variability in the association between behavioral time series. Psychological methods, 7(3), 338.

Moulder, R., Boker, S., Ramseyer, F., & Tschacher, W. (2018). Determining synchrony between behavioral time series: An application of surrogate data generation for establishing falsifiable null-hypotheses. Psychological Methods. 23:4 pp 757–773
}

\seealso{
\code{\link{wccAggregate}} for the definition of the columns of the returned dataframe. 
}

\examples{
\donttest{
#Create two arrays of timeseries with dyadic dependence
array1 <- matrix(NA, nrow=10, ncol=500)
array2 <- matrix(NA, nrow=10, ncol=500)
for(i in 1:10) {
    array1[i,] <- sin(c(1:500)/runif(1, min=5, max=20)) + rnorm(500, mean=0, sd=.5)
    array2[i,] <- array1[i,] + rnorm(500, mean=0, sd=.5)
}
wccSurrogateDyads(inArray1=array1, inArray2=array2, wMax=50, tMax=50, wInc=1, tInc=1,
    Lsize=8, pspan=.25, type="Max", nSurrogates=20, method="c")

}
}
