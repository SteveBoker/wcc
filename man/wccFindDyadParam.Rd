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

\name{wccFindDyadParam}
\alias{wccFindDyadParam}

\title{wccFindDyadParam}

\description{
This function calculates the difference between surrogate and actual dyads for a range of possible arguments to \code{wccCalc} and \code{wccPeakPick}.
}

\usage{
wccFindDyadParam(inArray1=NA, inArray2=NA, wMaxvector=c(50), tMaxvector=c(50),
    wIncvector=c(1), tIncvector=c(1), Lsizevector=c(8), pspanvector=c(.25), type="Max",
    nSurrogates=NA, samplespersecond=1,
    method=c("c", "cumr", "cumc", "cumcuda", "r"), embedD=9, ...)
}

\arguments{
\item{inArray1}{An N by T matrix of numeric values representing N timeseries of length T. These are the first set of time series which are randomly paired with timeseries from \sQuote{inArray2} and on which windowed cross correlation are calculated.  Must be of the same order N by T as \sQuote{inArray2}. Note that the only difference between \sQuote{inArray1} and \sQuote{inArray2} is which is treated as positive lag and which is negative lag in \code{wccPlot}.  When \sQuote{inArray1} is leading \sQuote{inArray2}, the peak correlation closest to zero in negative in \code{wccPlot}.  }
\item{inArray2}{An N by T matrix of numeric values representing N timeseries of length T. These are the second set of time series which are randomly paired with timeseries from \sQuote{inArray1} and on which windowed cross correlation are calculated.  Must be of the same order N by T as \sQuote{inArray1}.   When \sQuote{inArray2} is leading \sQuote{inArray1}, the peak correlation closest to zero in positive in \code{wccPlot}.}
\item{wMaxvector}{A numeric vector indicating the number of samples in a window for each test (default=c(50)).  This must be greater or equal to 5 samples.  It is recommended that wMax be greater than 15.}
\item{tMaxvector}{A numeric vector indicating the maximum number of samples to lag the windows for each test  (default=c(50)).  The cross correlation is calculated for windows that are time lagged against each other, allowing for processes than require some time to evolve.  For instance, if the two time series are from motion capture of two individuals conversing with one another, the speaker's movements are likely to occur prior to the responses of the listener.  tMax should be set to be greater than the number of samples that represents the maximum time lag that is likely to occur. }
\item{wIncvector}{A numeric vector indicating the number of samples incremented between successive windows for each test  (default=c(1)).  If the sample rate is very high and the time series is long, wInc may be increased.  However, \code{wccCalc}'s maximum sensitivity to nonstationary change in the time series occurs when wInc is 1.}
\item{tIncvector}{A numeric vector indicating the number of samples incremented between successive lags for each test  (default=c(1)).   If the sample rate is high and the maximum lag is high then tInc may be increased.  However, \code{wccCalc}'s maximum sensitivity to changes in lags occurs when the tInc is set to 1.}
\item{Lsizevector}{A numeric vector indicating  the width required for a peak to be identified by \code{wccPeakPick} for each test  (default=c(8)).  An Lsize that is too small may find spurious peaks that are just noise, whereas an Lsize that is too large may miss actual peaks.}
\item{pspanvector}{A numeric vector indicating the amount of smoothing applied by the loess function within \code{wccPeakPick} for each test  (default=c(.25)).}
\item{type}{A character string vector indicating whether a maximum, \sQuote{Max}, or minimum, \sQuote{Min}, should be detected by \code{wccPeakPick} for each test  (default="Max").}
\item{nSurrogates}{A numeric indicating the number of surrogates to generate (default=NA).  If there are N timeseries in \sQuote{inArray1}, then the maximum number of unique surrogates that can be generated is N*(N-1)/2}
\item{samplespersecond}{A numeric indicating the sampling rate of the time series (default=1).  This will determine the units of the first and second derivative aggregate statistics.}
\item{method}{Backend passed through to \code{\link{wccCalc}}: \code{"c"} (default, original compiled windcross), \code{"r"} (original pure-R loop), \code{"cumc"} (cumulative-sum C backend, which activates the batched code path), \code{"cumr"} (cumulative-sum pure-R reference), or \code{"cumcuda"} (cumulative-sum CUDA backend, also batched; requires a CUDA-enabled build).  See \code{\link{wccCalc}} for details and the missing-data caveat for the \code{cum*} methods.}
\item{embedD}{A numeric indicating the embedding dimension passed to \code{wccGLLAEmbed} to calculate derivatives of the peak picking lag timeseries.  The default (9) will provide a medium smooth at the cost of some loss of resolution in the timing of jumps in the second derivative.  The smallest legal value is 5, which provides less smooth and the best time resolution.  If your sampling rate is low, then a smaller value of \code{embedD} is useful, while if your sampling rate is high, a larger value of \code{embedD} will provide better noise rejection by smoothing.}
\item{...}{Reserved for the deprecated \code{windcross=TRUE/FALSE} argument, mapped to \code{method="c"} / \code{method="r"} with a warning.}
}

\details{
This function iterates over vectors of possible argument values given to \code{wccCalc} and \code{wccPeakPick} and for each combination of arguments returns the difference between the surrogate dyadic null hypothesis distribution and the true distribution of each of the statistics returned by wccAggregate.  The intended use of this function is to explore the space of possible argument values using a small hold-out set of dyadic timeseries in order to inform the selection of parameters to use for the full dataset.  It is recommended that between 10 and 20 dyadic timeseries be used as the hold-out set to input to this function.  The two arguments \sQuote{wMax} and \sQuote{tMax} are the most commonly explored parameters since they tend to have the most impact on results.  

Results from combinations of parameters can be visualized using the wccHeatmap function.  For any particular hypothesis, one or more aggregate statistics may be most useful in assessing the association between a dyad's time series.  Thus, one may wish to focus in on that statistic when choosing parameters.  wccHeatmapwMaxtMax allows the visualization of the effect of selecting the  \sQuote{wMax} and \sQuote{tMax} parameters on the statistic of interest.

Note that surrogate and real dyads will be calculated for all combinations of parameters in the argument vectors.  Thus, if one wishes to explore the effects of 8 choices of  \sQuote{wMax} and  8 choices of \sQuote{tMax}, a total of 64 surrogate analyses will be run.  If one uses a holdout set of 20 and asks for 100 surrogates, this will result in 7,680 calls to wccAggregate.  Such an exploration of the parameter space can take hours to run.  Add to that 2 choices of \sQuote{wInc} and  2 choices of \sQuote{tInc} and you may be talking about days of processing.  Be careful when choosing how many combinations to explore at one time.
}
    
\value{
Returns a dataframe with some of the columns from the wccAggregate function, some columns with redefined meanings, and some new columns.  

For each aggregated statistic a column giving the proportion of real dyads' statistics are greater than .95 or less than .05 of the values in the surrogate statistics distribution.  For each aggregated statistic a column giving the two sided Kolmogorov-Smirnov Test probabilities comparing the surrogate and real dyad distributions.  These probabilities are those returned by a call to ks.test() with alternative="two.sided".  For each aggregated statistic a column giving the quantile difference scores for the difference between the surrogate distribution and the real distribution of dyadic timeseries.  These quantile differences are the mean of the difference between calls to quantile() with probs=seq(.1,.9, by=.1) for the real dyads and the surrogate dyads.  Each row in the dataframe presents the results from one combination of chosen parameters.

The returned data frame is defined as follows.  
\tabular{ll}{
\sQuote{wMax}   \tab  The value of wMax in the argument passed to \code{wccCalc}.\cr

\sQuote{tMax}   \tab  The value of tMax in the argument passed to \code{wccCalc}. \cr

\sQuote{wInc}   \tab  The value of wInc in the argument passed to \code{wccCalc}. \cr

\sQuote{tInc}   \tab   The value of tInc in the argument passed to \code{wccCalc}.\cr

\sQuote{Lsize}   \tab  The value of Lsize in the argument passed to \code{wccPeakPick}. \cr

\sQuote{pspan}   \tab  The value of pspan in the argument passed to \code{wccPeakPick}. \cr

\sQuote{type}   \tab   The value of type as a numeric translation of the type argument passed to \code{wccPeakPick}.  This is +1 when type = \sQuote{Max} and -1 when type = \sQuote{Min}.\cr

\sQuote{samples}   \tab  Length of the time series inSeries1 and inSeries2\cr

\sQuote{windows}   \tab  The number of windows returned by \code{wccCalc}.\cr

\sQuote{pctMissing}   \tab  Proportion of missing values in the matrix returned by \code{wccCalc}.\cr

\sQuote{pctMissingWindows}   \tab  Proportion of windows (i.e., rows) in the matrix returned by \code{wccCalc} that include at least one missing value.\cr

\sQuote{maxMean}   \tab The proportion of the real data that are either greater than 95\%  less than 5\% of the surrogate distribution of means of the Fischer's Z transform of the maximum correlation vector returned returned by \code{wccPeakPick}.  Note that if type=\sQuote{Min} then this value is the Fischer's Z transform of the minimum correlation vector. \cr

\sQuote{maxVar}   \tab   The proportion of the real data that are either greater than 95\%  less than 5\% of the surrogate distribution of variance of the Fischer's Z transform of the maximum correlation vector returned returned by \code{wccPeakPick}.  Note that if type=\sQuote{Min} then this value is the Fischer's Z transform of the minimum correlation vector. \cr

\sQuote{totalMean}   \tab  The proportion of the real data that are either greater than 95\%  less than 5\% of the surrogate distribution of mean of the Fischer's Z transform of all values in the matrix returned by \code{wccCalc}.\cr

\sQuote{totalVar}   \tab  The proportion of the real data that are either greater than 95\%  less than 5\% of the surrogate distribution of variance of the Fischer's Z transform of all values in the matrix returned by \code{wccCalc} .\cr

\sQuote{zeroLagMean}   \tab  The proportion of the real data that are either greater than 95\%  less than 5\% of the surrogate distribution of mean of the Fischer's Z transform of all values in the zero lag column of the matrix returned by \code{wccCalc}.\cr

\sQuote{zeroLagVar}   \tab The proportion of the real data that are either greater than 95\%  less than 5\% of the surrogate distribution of variance of the Fischer's Z transform of all values in the zero lag column of the matrix returned by \code{wccCalc}. \cr

\sQuote{lagMean}   \tab  The proportion of the real data that are either greater than 95\%  less than 5\% of the surrogate distribution of mean of the lag vector returned by \code{wccPeakPick}.\cr

\sQuote{lagVar}   \tab  The proportion of the real data that are either greater than 95\%  less than 5\% of the surrogate distribution of variance of the lag vector returned by \code{wccPeakPick}.\cr

\sQuote{dlagMean}   \tab The proportion of the real data that are either greater than 95\%  less than 5\% of the surrogate distribution of mean of the first derivative of the lag vector returned by \code{wccPeakPick}. \cr

\sQuote{dlagVar}   \tab  The proportion of the real data that are either greater than 95\%  less than 5\% of the surrogate distribution of variance of the first derivative of the lag vector returned by \code{wccPeakPick}.\cr

\sQuote{d2lagMean}   \tab The proportion of the real data that are either greater than 95\%  less than 5\% of the surrogate distribution of mean of the second derivative of the lag vector returned by \code{wccPeakPick}. \cr

\sQuote{d2lagVar}   \tab  The proportion of the real data that are either greater than 95\%  less than 5\% of the surrogate distribution of variance of the second derivative of the lag vector returned by \code{wccPeakPick}.\cr

\sQuote{maxMean}   \tab  Mean of the Fischer's Z transform of the maximum correlation vector returned returned by \code{wccPeakPick}.  Note that if type=\sQuote{Min} then this value is the Fischer's Z transform of the minimum correlation vector. \cr

\sQuote{maxVar}   \tab   Variance of the Fischer's Z transform of the maximum correlation vector returned returned by \code{wccPeakPick}.  Note that if type=\sQuote{Min} then this value is the Fischer's Z transform of the minimum correlation vector. \cr

\sQuote{totalMeanKS}   \tab  The probability of the two tailed Kolmogorov-Smirnov test for the difference between the surrogate distribution and real dyad distribution of the mean of the Fischer's Z transform of all values in the matrix returned by \code{wccCalc}.\cr

\sQuote{totalVarKS}   \tab  The probability of the two tailed Kolmogorov-Smirnov test for the difference between the surrogate distribution and real dyad distribution of the variance of the Fischer's Z transform of all values in the matrix returned by \code{wccCalc} .\cr

\sQuote{zeroLagMean}   \tab  The probability of the two tailed Kolmogorov-Smirnov test for the difference between the surrogate distribution and real dyad distribution of the mean of the Fischer's Z transform of all values in the zero lag column of the matrix returned by \code{wccCalc}.\cr

\sQuote{zeroLagVarKS}   \tab The probability of the two tailed Kolmogorov-Smirnov test for the difference between the surrogate distribution and real dyad distribution of the variance of the Fischer's Z transform of all values in the zero lag column of the matrix returned by \code{wccCalc}. \cr

\sQuote{lagMeanKS}   \tab  The probability of the two tailed Kolmogorov-Smirnov test for the difference between the surrogate distribution and real dyad distribution of the mean of the lag vector returned by \code{wccPeakPick}.\cr

\sQuote{lagVarKS}   \tab  The probability of the two tailed Kolmogorov-Smirnov test for the difference between the surrogate distribution and real dyad distribution of the variance of the lag vector returned by \code{wccPeakPick}.\cr

\sQuote{dlagMeanKS}   \tab The probability of the two tailed Kolmogorov-Smirnov test for the difference between the surrogate distribution and real dyad distribution of the mean of the first derivative of the lag vector returned by \code{wccPeakPick}. \cr

\sQuote{dlagVarKS}   \tab  The probability of the two tailed Kolmogorov-Smirnov test for the difference between the surrogate distribution and real dyad distribution of the variance of the first derivative of the lag vector returned by \code{wccPeakPick}.\cr

\sQuote{d2lagMeanKS}   \tab The probability of the two tailed Kolmogorov-Smirnov test for the difference between the surrogate distribution and real dyad distribution of the mean of the second derivative of the lag vector returned by \code{wccPeakPick}. \cr

\sQuote{d2lagVarKS}   \tab  The probability of the two tailed Kolmogorov-Smirnov test for the difference between the surrogate distribution and real dyad distribution of the variance of the second derivative of the lag vector returned by \code{wccPeakPick}.\cr

\sQuote{maxMeanQdiff}   \tab  The quantile difference (using deciles) between the surrogate and the real distribution of the mean of the Fischer's Z transform of the maximum correlation vector returned returned by \code{wccPeakPick}.  Note that if type=\sQuote{Min} then this value is the Fischer's Z transform of the minimum correlation vector. \cr

\sQuote{maxVarQdiff}   \tab   The quantile difference (using deciles) between the surrogate and the real distribution of the Variance of the Fischer's Z transform of the maximum correlation vector returned returned by \code{wccPeakPick}.  Note that if type=\sQuote{Min} then this value is the Fischer's Z transform of the minimum correlation vector. \cr

\sQuote{totalMeanQdiff}   \tab  The quantile difference (using deciles) between the surrogate and the real distribution of the mean of the Fischer's Z transform of all values in the matrix returned by \code{wccCalc}.\cr

\sQuote{totalVarQdiff}   \tab  The quantile difference (using deciles) between the surrogate and the real distribution of the variance of the Fischer's Z transform of all values in the matrix returned by \code{wccCalc} .\cr

\sQuote{zeroLagMeanQdiff}   \tab  The quantile difference (using deciles) between the surrogate and the real distribution of the mean of the Fischer's Z transform of all values in the zero lag column of the matrix returned by \code{wccCalc}.\cr

\sQuote{zeroLagVarQdiff}   \tab The quantile difference (using deciles) between the surrogate and the real distribution of the variance of the Fischer's Z transform of all values in the zero lag column of the matrix returned by \code{wccCalc}. \cr

\sQuote{lagMeanQdiff}   \tab  The quantile difference (using deciles) between the surrogate and the real distribution of the mean of the lag vector returned by \code{wccPeakPick}.\cr

\sQuote{lagVarQdiff}   \tab  The quantile difference (using deciles) between the surrogate and the real distribution of the variance of the lag vector returned by \code{wccPeakPick}.\cr

\sQuote{dlagMeanQdiff}   \tab The quantile difference (using deciles) between the surrogate and the real distribution of the mean of the first derivative of the lag vector returned by \code{wccPeakPick}. \cr

\sQuote{dlagVarQdiff}   \tab  The quantile difference (using deciles) between the surrogate and the real distribution of the variance of the first derivative of the lag vector returned by \code{wccPeakPick}.\cr

\sQuote{d2lagMeanQdiff}   \tab The quantile difference (using deciles) between the surrogate and the real distribution of the mean of the second derivative of the lag vector returned by \code{wccPeakPick}. \cr

\sQuote{d2lagVarQdiff}   \tab  The quantile difference (using deciles) between the surrogate and the real distribution of the variance of the second derivative of the lag vector returned by \code{wccPeakPick}.\cr
}

}

\references{
Boker, S. M., Rotondo, J. L., Xu, M., & King, K. (2002). Windowed cross-correlation and peak picking for the analysis of variability in the association between behavioral time series. Psychological methods, 7(3), 338.

William J. Conover (1971).  Practical Nonparametric Statistics. New York: John Wiley & Sons.  Pages 295-301 (one-sample Kolmogorov test), 309-314 (two-sample Smirnov test).
}

\seealso{
\code{\link{wccHeatmap}} visualizes results selected aggregated statistics for combinations of parameters.  See \code{\link{wccAggregate}} for definitions of the statistics that are compared between surrogate and real dyads.  See \code{\link{ks.test}} for definition of the two tailed Kolmogorov-Smirnov Test.
}

\examples{
\donttest{
# Create two arrays of timeseries with dyadic dependence
array1 <- matrix(NA, nrow=10, ncol=500)
array2 <- matrix(NA, nrow=10, ncol=500)
for(i in 1:10) {
    array1[i,] <- sin(c(1:500)/runif(1, min=5, max=20)) + rnorm(500, mean=0, sd=.5)
    array2[i,] <- array1[i,] + rnorm(500, mean=0, sd=.5)
}
# Select parameters to explore
wMaxVec <- c(50,100)
tMaxVec <- c(25,50)
wccFindDyadParam(inArray1=array1, inArray2=array2, wMaxvector=wMaxVec, tMaxvector=tMaxVec,
    nSurrogates=30, samplespersecond=1, method="c")
}
}
