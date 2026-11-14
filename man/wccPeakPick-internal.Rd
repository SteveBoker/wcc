\name{wccPeakPick-internal}
\alias{wccSmoothMatrix}
\alias{wccPeakSearch}

\title{Internal helpers for wccPeakPick}

\description{
Internal helpers used by \code{\link{wccPeakPick}} and \code{\link{wccPeakPickBatch}}.  \code{wccSmoothMatrix} builds and caches the loess-equivalent smoother matrix for a given column length and \code{pspan}; \code{wccPeakSearch} scans a smoothed row for a peak given a bounding window.  Not intended for direct use.
}

\keyword{internal}
