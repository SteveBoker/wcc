# ---------------------------------------------------------------------
# Program: testWccSymmetry.R
#
# Tests for wccZeroChangeDirections and wccSymmetryEpochs.
#
# 1. Constructed matrices: a rigidly translating lag profile must have
#    near-zero dispersion and recover its drift rate; stationary and
#    drifting epochs must be split when driftTolerance is set.
# 2. wccCalc on simulated series with a known coupled epoch (fixed lag)
#    and a known drifting-lag epoch: detected formation and breaking
#    times must fall near the true boundaries, and uncoupled stretches
#    must not be marked symmetric.
# ---------------------------------------------------------------------

library(wcc)

set.seed(261106)

# ---- 1. Constructed matrices ----------------------------------------

lags <- -20:20
nT <- 200
bump <- function(center) .8 * exp(-(lags - center)^2 / (2 * 3^2))
k <- dnorm(-3:3, 0, 1.5); k <- k / sum(k)
noise <- matrix(rnorm(nT * length(lags)), nT)
noise <- t(apply(noise, 1, function(r) stats::filter(r, k, circular=TRUE)))
noise <- apply(noise, 2, function(cl) stats::filter(cl, k, circular=TRUE))
noise <- noise / sd(noise) * .3

m <- noise
for (i in 51:100)  m[i, ] <- bump(5) + .03 * noise[i, ]
for (i in 101:150) m[i, ] <- bump(-10 + .4 * (i - 101)) + .03 * noise[i, ]

z <- wccZeroChangeDirections(m, method="gradient")
stopifnot(median(z$dispersion[55:95])   < .05)
stopifnot(median(z$dispersion[105:145]) < .05)
stopifnot(median(z$dispersion[c(5:45, 155:195)]) > .4)
stopifnot(abs(median(z$driftRate[55:95])) < .02)
stopifnot(abs(median(z$driftRate[105:145]) - .4) < .02)

e <- wccSymmetryEpochs(z)$epochs
stopifnot(nrow(e) == 1, abs(e$formationWindow - 51) <= 2, abs(e$breakingWindow - 151) <= 2)

e <- wccSymmetryEpochs(z, driftTolerance=.1)$epochs
stopifnot(nrow(e) == 2)
stopifnot(abs(e$formationWindow - c(51, 101)) <= 2)
stopifnot(abs(e$breakingWindow  - c(101, 151)) <= 2)
stopifnot(abs(e$medianDriftRate - c(0, .4)) < .02)

# drift rate is reported in lag seconds per second
z2 <- wccZeroChangeDirections(m, method="gradient", wInc=10, tInc=2, samplespersecond=5)
stopifnot(abs(median(z2$driftRate[105:145]) - .4 * 2 / 10) < .005)

# shift method: same epochs, drift recovered from the best uniform shift
zs <- wccZeroChangeDirections(m, timeScale=5, maxShift=8)
stopifnot(median(zs$dispersion[60:90])   < .05)
stopifnot(median(zs$dispersion[110:140]) < .05)
stopifnot(median(zs$dispersion[c(10:40, 160:190)]) > .3)
stopifnot(abs(median(zs$driftRate[60:90])) < .02)
stopifnot(abs(median(zs$driftRate[110:140]) - .4) < .02)
stopifnot(all(is.na(zs$dispersion[c(1:5, 196:200)])))

# ---- 2. wccCalc on simulated series ---------------------------------

N <- 3000
x <- as.numeric(arima.sim(list(ar=.9), N))
y <- as.numeric(arima.sim(list(ar=.9), N))
# samples 800-1500: y follows x at lag 6
for (i in 800:1500) y[i] <- x[i - 6] + .3 * rnorm(1)
# samples 1800-2400: lag drifts from 0 to 15
for (i in 1800:2400) y[i] <- x[i - round((i - 1800) / 40)] + .3 * rnorm(1)

wMax <- 100; tMax <- 20; wInc <- 10
W <- wccCalc(x, y, wMax=wMax, tMax=tMax, wInc=wInc, tInc=1)
z <- wccZeroChangeDirections(W, wMax=wMax, tMax=tMax, wInc=wInc, method="gradient", timeSmooth=wMax/wInc)
e <- wccSymmetryEpochs(z, minDuration=10, maxGap=5)$epochs
print(e)

# Window-level separation: coupled windows (away from the boundaries)
# have low dispersion, uncoupled windows have high dispersion.
coupled  <- (z$elapsedTime >= 900 & z$elapsedTime <= 1450) |
            (z$elapsedTime >= 1900 & z$elapsedTime <= 2350)
uncoupled <- z$elapsedTime < 700 | z$elapsedTime > 2550
stopifnot(median(z$dispersion[coupled]) < .2)
stopifnot(median(z$dispersion[uncoupled]) > .3)

# Both coupled epochs begin near their true start.
near <- function(a, b, tol) any(abs(a - b) <= tol, na.rm=TRUE)
stopifnot(near(e$formationTime, 800,  1.5 * wMax))
stopifnot(near(e$formationTime, 1800, 1.5 * wMax))

# ---- 3. Clean sines with noise intervals ------------------------------
# Two sines of slightly different frequency drift uniformly in relative
# phase (symmetric); replacing series 2 with noise breaks symmetry for
# windows whose span (wMax + tMax samples) reaches into the noise.

set.seed(1)
tSeries1 <- sin(c(1:1000)/10) + rnorm(1000, mean=0, sd=.005)
tSeries2 <- sin(c(11:1010)/11) + rnorm(1000, mean=0, sd=.005)
tSeries2[201:300] <- rnorm(100, mean=0, sd=1)
tSeries2[501:600] <- rnorm(100, mean=0, sd=1)
W <- wccCalc(tSeries1, tSeries2, wMax=100, tMax=50, wInc=1, tInc=1)
sym <- wccSymmetryEpochs(W, wMax=100, tMax=50, method="gradient", minDuration=20, maxGap=20)
e <- sym$epochs
print(e)
# symmetry breaks as windows reach the noise at 201 and 501 ...
stopifnot(near(e$breakingTime, 201, 25), near(e$breakingTime, 501, 25))
# ... and re-forms once they clear it (about 300 + 150 and 600 + 150)
stopifnot(near(e$formationTime, 450, 25), near(e$formationTime, 750, 25))
# no epoch covers the cores of the disrupted intervals
w <- sym$windows
stopifnot(!any(w$symmetric[w$elapsedTime >= 230 & w$elapsedTime <= 400]))
# the uniform phase drift is recovered: (1/10 - 1/11) / (1/10) = 0.091 lag/sample
stopifnot(abs(abs(median(e$medianDriftRate)) - 1/11) < .02)

# shift method with a data-driven threshold finds the same breaks
W <- wccCalc(tSeries1, tSeries2, wMax=30, tMax=50, wInc=1, tInc=1)
sym <- wccSymmetryEpochs(W, wMax=30, tMax=50, threshold="auto")
e <- sym$epochs
print(e[, c("startTime", "endTime", "medianDriftRate")])
stopifnot(near(e$endTime, 200, 30), near(e$endTime, 500, 30))
stopifnot(near(e$startTime, 300, 40), near(e$startTime, 600, 40))
w <- sym$windows
stopifnot(!any(w$symmetric[w$elapsedTime >= 250 & w$elapsedTime <= 300]))

# ---- 4. Shift method defaults on a noisy dyad with known epochs -------
# (same construction as the generator SymmetryEpochTest-261111.R)

set.seed(3)
n <- 1000
x <- rnorm(n)
y <- rnorm(n)
s <- sin(1:n/10) + rnorm(n, sd=.5)
x[101:250] <- s[101:250];  y[101:250] <- s[(101:250) + 8]
s1 <- sin(1:n/9) + rnorm(n, sd=.5);  s2 <- sin(1:n/11) + rnorm(n, sd=.5)
x[301:480] <- s1[301:480]; y[301:480] <- s2[(301:480) + 5]
s3 <- sin(1:n/10) + rnorm(n, sd=.5); r3 <- .5 * s3 + rnorm(n, sd=.5)
x[531:700] <- s3[531:700]; y[531:700] <- r3[(531:700) - 6]
sym <- wccSymmetryEpochs(wccCalc(x, y, wMax=30, tMax=30), wMax=30, tMax=30)
e <- sym$epochs
print(e[, c("startTime", "endTime", "peakCorrelation", "peakLagStability", "peakSharpness")])
# the lagged copy and the drifting sines are found within 20 samples
stopifnot(abs(e$startTime[1:2] - c(101, 301)) <= 20)
stopifnot(abs(e$endTime[1:2] - c(250, 480)) <= 20)
# the weak regression is detected over most of its extent (it may split)
covered <- sum(sapply(531:700, function(t) any(e$startTime <= t & e$endTime >= t)))
stopifnot(covered / 170 > .6)
# nothing is detected in the noise after sample 720
stopifnot(!any(e$endTime > 720 & e$startTime > 720))
# lagged copy: near-perfect, stable, sharp peak at lag 8
stopifnot(e$peakCorrelation[1] > .95, e$peakLagStability[1] > .9, abs(abs(e$peakLag[1]) - 8) <= 1)
# different periods: the peak moves
stopifnot(e$peakLagStability[2] < .6)
# epochs do not overlap
stopifnot(all(e$startTime[-1] > e$endTime[-nrow(e)]))

cat("testWccSymmetry: all tests passed\n")
