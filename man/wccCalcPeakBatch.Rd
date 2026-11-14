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

\name{wccCalcPeakBatch}
\alias{wccCalcPeakBatch}

\title{Fused Batched WCC and Peak Pick}

\description{
Computes the Windowed Cross Correlation grids for a list of series pairs and immediately peak-picks each grid, in one call.  Equivalent to \code{wccCalcBatch} followed by an \code{NA -> 0} substitution and \code{wccPeakPickBatch}, but with \code{method="cumcuda"} the intermediate grids never leave the GPU: only the small per-pair index and value matrices are transferred back to the host.  For GPU workloads this avoids the dominant cost of the two-step pipeline (a device-host-device round trip of the full grid array).
}

\usage{
wccCalcPeakBatch(seriesArray1, seriesArray2, pairs=NULL,
                 wMax=50, tMax=50, wInc=1, tInc=1,
                 Lsize=8, pspan=.25, type="Max",
                 method=c("cumcuda", "cumc"),
                 precision=c("double", "single"))
}

\arguments{
\item{seriesArray1}{A \code{D1 x n} numeric matrix; each row is one series.}
\item{seriesArray2}{A \code{D2 x n} numeric matrix; each row is one series.  Must have the same number of columns as \code{seriesArray1}.}
\item{pairs}{A \code{P x 2} integer matrix; row \code{p} selects the pair \code{(seriesArray1[pairs[p,1],], seriesArray2[pairs[p,2],])}.  When \code{NULL} (default), it is set to \code{cbind(1:D, 1:D)}, i.e. the real dyads.}
\item{wMax}{Window size in samples (default=50).  Must be at least 5.}
\item{tMax}{Maximum lag in samples (default=50).}
\item{wInc}{Window increment in samples (default=1).}
\item{tInc}{Lag increment in samples (default=1).}
\item{Lsize}{Look-ahead size of the peak search (default=8); see \code{\link{wccPeakPick}}.}
\item{pspan}{Span of the loess smoother (default=.25); see \code{\link{wccPeakPick}}.}
\item{type}{\code{"Max"} to find peaks (default) or \code{"Min"} to find valleys.}
\item{method}{\code{"cumcuda"} (default) runs both stages on the GPU with the grid resident on the device (requires a CUDA-enabled build); \code{"cumc"} composes the two-step CPU pipeline.}
\item{precision}{Floating-point precision of the GPU computation; only valid with \code{method="cumcuda"}.  \code{"double"} (default) or \code{"single"}.}
}

\details{
Zero-variance windows (which \code{wccCalcBatch} reports as \code{NA}) are treated as correlation 0 before smoothing, matching the conventional \code{g[is.na(g)] <- 0} preprocessing of the two-step pipeline.  Missing data in the input series is not supported and raises an error.
}

\value{
A list with elements \code{index} and \code{value}, each an \code{nRow x P} matrix: the peak lag position and smoothed peak value per window row and pair, \code{NA} where the peak fails the boundary criterion.
}

\seealso{
\code{\link{wccCalcBatch}}, \code{\link{wccPeakPickBatch}}
}

\examples{
D <- 3
n <- 500
arr1 <- t(sapply(1:D, function(i) sin(c(1:n)/10) + rnorm(n, sd=.5)))
arr2 <- t(sapply(1:D, function(i) sin(c(1:n)/12) + rnorm(n, sd=.5)))
pp <- wccCalcPeakBatch(arr1, arr2, wMax=50, tMax=50, method="cumc")
dim(pp$index)
}
