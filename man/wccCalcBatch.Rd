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

\name{wccCalcBatch}
\alias{wccCalcBatch}

\title{Batched Windowed Cross Correlation}

\description{
Computes the Windowed Cross Correlation grid for a list of series pairs in one call, reusing per-series prefix sums across pairs.  Real dyads and surrogate pairings drawn from the same series pool can therefore be evaluated together cheaply.
}

\usage{
wccCalcBatch(seriesArray1, seriesArray2, pairs=NULL,
             wMax=50, tMax=50, wInc=1, tInc=1, method=c("cumc", "cumcuda"),
             precision=c("double", "single"))
}

\arguments{
\item{seriesArray1}{A \code{D1 x n} numeric matrix; each row is one series.  When \code{pairs} is not supplied, \code{seriesArray1} and \code{seriesArray2} must have the same number of rows.}
\item{seriesArray2}{A \code{D2 x n} numeric matrix; each row is one series.  Must have the same number of columns as \code{seriesArray1}.}
\item{pairs}{A \code{P x 2} integer matrix; row \code{p} selects the pair \code{(seriesArray1[pairs[p,1],], seriesArray2[pairs[p,2],])}.  When \code{NULL} (default), it is set to \code{cbind(1:D, 1:D)}, i.e. the real dyads.}
\item{wMax}{Window size in samples (default=50).  Must be at least 5.}
\item{tMax}{Maximum lag in samples (default=50).}
\item{wInc}{Window increment in samples (default=1).}
\item{tInc}{Lag increment in samples (default=1).}
\item{method}{Backend for the batched computation: \code{"cumc"} (default, batched cumulative-sum C backend) or \code{"cumcuda"} (batched CUDA backend, requires a CUDA-enabled build).}
\item{precision}{Floating-point precision of the GPU computation; only valid with \code{method="cumcuda"}.  \code{"double"} (default) computes in FP64.  \code{"single"} computes in FP32, which is substantially faster on consumer GPUs at the cost of correlations accurate to roughly 1e-5 rather than machine precision.  Results are always returned as R doubles.}
}

\details{
Missing data is not supported: any \code{NA} in \code{seriesArray1} or \code{seriesArray2} raises an error.

The C backend computes per-series prefix sums once and reuses them across all pairs, so the batched call scales with the number of unique series rather than the number of pairs.  This makes the function particularly efficient when many surrogate pairings are drawn from a small pool of series, as in \code{\link{wccSurrogateDyads}} and \code{\link{wccFindDyadParam}}.
}

\value{
An \code{nRow x nCol x P} numeric array.  Slab \code{[,,p]} equals \code{wccCalc(seriesArray1[pairs[p,1],], seriesArray2[pairs[p,2],], wMax=wMax, tMax=tMax, wInc=wInc, tInc=tInc, method="cumc")}.
}

\seealso{
\code{\link{wccCalc}} for the single-dyad interface and the definition of the returned grid.
}

\examples{
D <- 3
n <- 500
arr1 <- t(sapply(1:D, function(i) sin(c(1:n)/10) + rnorm(n, sd=.5)))
arr2 <- t(sapply(1:D, function(i) sin(c(1:n)/12) + rnorm(n, sd=.5)))
grids <- wccCalcBatch(arr1, arr2, wMax=50, tMax=50)
dim(grids)
}
