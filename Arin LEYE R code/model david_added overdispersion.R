model
{
  for (i in 1:n) {
    count[i] ~ dpois(mu[i])
    mu[i] <- exp(bw[w[i]] + by[y[i]] + bet * et[i] + bed * ed[i] + e[i])
    e[i] ~ dnorm(0, tau.e)
    countsim[i] ~ dpois(mu[i])
  }

  for (i in 1:nw) {
    bw[i] <- bwc *
      wc[i] +
      bwa * wa[i] +
      bwt[wt[i]] +
      bwh[wh[i]] +
      bwu[wu[i]] +
      ew[i]
    ew[i] ~ dnorm(0, tau.ew)
  }

  for (i in 1:ny) {
    by[i] <- b0 + b1 * i + ey[i]
    ey[i] ~ dnorm(0, tau.ey)
  }

  bet ~ dnorm(0, 0.001)
  bed ~ dnorm(0, 0.001)

  bwc ~ dnorm(0, 0.001)
  bwa ~ dnorm(0, 0.001)

  bwt[1] <- 0
  bwh[1] <- 0
  bwu[1] <- 0

  for (i in 2:nwt) {
    bwt[i] ~ dnorm(0, 0.001)
  }
  for (i in 2:nwh) {
    bwh[i] ~ dnorm(0, 0.001)
  }
  for (i in 2:nwu) {
    bwu[i] ~ dnorm(0, 0.001)
  }

  b0 ~ dnorm(0, 0.001)
  b1 ~ dnorm(0, 0.001)

  tau.e <- 1 / (sd.e * sd.e)
  sd.e ~ dunif(0, 10)

  tau.ew <- 1 / (sd.ew * sd.ew)
  sd.ew ~ dunif(0, 10)

  tau.ey <- 1 / (sd.e * sd.ey)
  sd.ey ~ dunif(0, 10)
}
