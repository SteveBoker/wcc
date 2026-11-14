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

\name{wccPeakPickBatch}
\alias{wccPeakPickBatch}

\title{Batched Peak Picking Across a Stack of WCC Grids}

\description{
Runs peak picking over a 3-D array of Windowed Cross Correlation grids and returns the picked index and value for every slab.  Intended for consuming the array returned by \code{\link{wccCalcBatch}}.
}

\usage{
wccPeakPickBatch(grids, Lsize=8, pspan=.25, type="Max", method=c("cumc", "cumcuda"),
                 precision=c("double", "single"))
}

\arguments{
\item{grids}{A 3-dimensional numeric array of shape \code{nRow x nCol x P} whose slabs \code{grids[,,p]} are WCC grids to peak-pick.}
\item{Lsize}{An integer specifying the width required for a peak to be identified (default=8).}
\item{pspan}{A numeric specifying the amount of smoothing applied by the shared matrix smoother (default=0.25).}
\item{type}{Either \code{"Max"} or \code{"Min"}, selecting peaks or valleys respectively.}
\item{method}{Backend selector.  \code{"cumc"} (default) runs the CPU path per slab, reusing the cached matrix smoother.  \code{"cumcuda"} smooths all slabs with one batched GEMM on the GPU and runs the search kernel there (requires a CUDA-enabled build); NA cells must be replaced before calling.}
\item{precision}{Floating-point precision of the GPU computation; only valid with \code{method="cumcuda"}.  \code{"double"} (default) computes in FP64.  \code{"single"} smooths and searches in FP32, which is substantially faster on consumer GPUs; peak values are then accurate to roughly 1e-5 and near-tied peak positions may occasionally shift by one half-step.  Results are always returned as R doubles.}
}

\value{
A list with components \code{index} and \code{value}, each an \code{nRow x P} matrix.  Column \code{p} of \code{index}/\code{value} contains the picked index/value for the \code{p}-th slab of \code{grids}, matching the maximum-side output of \code{\link{wccPeakPick}} on the same slab.
}

\seealso{
\code{\link{wccPeakPick}} for the single-slab interface and details of the smoother; \code{\link{wccCalcBatch}} for the batched WCC producer whose output feeds this function.
}

\examples{
n <- 500
g1 <- wccCalc(sin(c(1:n)/10) + rnorm(n, sd=.5),
              sin(c(1:n)/12) + rnorm(n, sd=.5), wMax=25, tMax=25, method="cumc")
g2 <- wccCalc(sin(c(1:n)/11) + rnorm(n, sd=.5),
              sin(c(1:n)/13) + rnorm(n, sd=.5), wMax=25, tMax=25, method="cumc")
g1[is.na(g1)] <- 0
g2[is.na(g2)] <- 0
grids <- array(NA_real_, dim=c(nrow(g1), ncol(g1), 2))
grids[,,1] <- g1
grids[,,2] <- g2
wccPeakPickBatch(grids, Lsize=8, pspan=.25, type="Max", method="cumc")
}
