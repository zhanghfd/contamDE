###  different estimation of deviations and degrees of freedom
ebayes <- function(var.res, dg,winsor.tail.p = c(0.05, 0.1), robust = TRUE) {

  n <- length(var.res)
  df1 <- dg
  prob <- winsor.tail.p
  prob[2] <- 1 - winsor.tail.p[2]
  xq <- quantile(var.res, prob)

  ## non robust  (without outliers)
  if(!robust) {
    zg <- log(var.res)  ## fisher's z distribution
    z.mean <- mean(zg)
    z.var <- var(zg)

    fund0 <- function(x) {
      trigamma(x / 2) + trigamma(dg / 2) - z.var
    }

    d0 <- uniroot(fund0, c(0.1,500))$root

    ## the estimation of s0^2

    log.s02 <- digamma(d0 / 2) - log(d0 / dg) - digamma(dg / 2) + z.mean
    s02 <- exp(log.s02)

    ## residual variance EB estimation of sigma_g^2
    df <- d0 + dg
    sigma.g2.hat=(d0 * s02 + dg * var.res) / df

    result = list(sigma.g2.hat = sigma.g2.hat, d0 = d0, s02 = s02)
    return(result)
  }

  ## robust (with outliers)
  win.s.g <- pmax(pmin(var.res, xq[2]), xq[1])
  zg <- log(win.s.g)
  z.mean <- mean(zg)
  z.var  <- var(zg)

  sub.fun <- function(x) x / (1 + x)
  sub.inv <- function(x) x / (1 - x)
  gauss.quad <- statmod::gauss.quad.prob(128, dist = "uniform")

  Win.F.Moments <- function(df1 = df1, df2 = df2) {
    fq <- qf(p = c(winsor.tail.p[1], 1 - winsor.tail.p[2]),
             df1 = df1, df2 = df2)
    zq <- log(fq)
    q <- sub.fun(fq)
    q21 <- q[2] - q[1]
    nodes <- q[1] + (q[2] - q[1]) * gauss.quad$nodes
    inv.nodes <- sub.inv(nodes)
    pdf.nodes <- df(inv.nodes, df1 = df1, df2 = df2)
    log.inv.nodes <- log(inv.nodes)
    h.mean <- log.inv.nodes * pdf.nodes / (1 - nodes) ^ 2
    mean <- q21 * sum(gauss.quad$weight * h.mean) + sum(zq * winsor.tail.p)
    h.var  <- (log.inv.nodes - mean) ^ 2 / (1 - nodes) ^ 2 * pdf.nodes
    var <- q21 * sum(gauss.quad$weight * h.var) +
      sum((zq - mean) ^ 2 * winsor.tail.p)

    list(mean = mean, var = var)
  }

  inf.mom <- Win.F.Moments(df1 = dg, df2 = Inf)
  ub.fun <- log(z.var/inf.mom$var)

  if(ub.fun <= 0) {
    df2 <- rep(Inf, n)
    log.s02 <- z.mean - inf.mom$mean
    s02 <- exp(log.s02)
    sigma.g2.hat <- rep(0, n)

    return(list(sigma.g2.hat = sigma.g2.hat, d0 = df2.corrected, s02 = s02))
  }

  non.robust <- ebayes(var.res, dg,winsor.tail.p = c(0.05, 0.1), robust = FALSE)
  fun <- function(x) {
    df2 <- sub.inv(x)
    mom <- Win.F.Moments(df1 = df1, df2 = df2)
    log(z.var/mom$var)
  }
  lbx <- sub.fun(non.robust$d0)
  lb.fun <- fun(lbx)

  if (lb.fun >= 0) {
    df2 <- non.robust$d0
  } else {
    u <- uniroot(fun, interval = c(lbx, 1), tol = 1e-08,f.lower = lb.fun, f.upper = ub.fun)
    df2 <- sub.inv(u$root)
  }

  mom <- Win.F.Moments(df1 = df1, df2 = df2)
  log.s02 <- z.mean - mom$mean
  s02 <- exp(log.s02)

  Fstat <- exp(log(var.res) - log.s02)

  log.p.value <- pf(Fstat, df1 = df1, df2 = df2, lower.tail = FALSE,log.p = TRUE)

  p.value <- exp(log.p.value)

  r <- rank(var.res)
  log.prior.p <- log(n - r + 0.5) - log(n)
  log.prob.not.outlier <- pmin(log.p.value - log.prior.p, 0)
  prob.not.outlier <- exp(log.prob.not.outlier)
  prob.outlier <- -expm1(log.prob.not.outlier)

  likelihood.fun <- function(x) {
    df(max(var.res)/s02, df1 = df1, df2 = x, log = T)
  }

  df2.outlier <- optimise(likelihood.fun, c(0, df2), maximum = TRUE)$maximum
  df2.corrected <- prob.not.outlier * df2 + prob.outlier * df2.outlier

    if (any(log.prob.not.outlier < 0)) {

      o <- order(log.p.value)
      df2.ordered <- df2.corrected[o]
      cummean.df2 <- cumsum(df2.ordered) / (1: n)
      min.index <- which.min(cummean.df2)

      df2.ordered[1:min.index] <- cummean.df2[min.index]
      df2.corrected[o] <- cummax(df2.ordered)
    } else {
      df2.outlier <- df2.outlier2 <- df2
      df2.corrected <- rep.int(df2, n)
    }

  df <- df2.corrected + dg
  sigma.g2.hat <- (df2.corrected * s02 + dg * var.res) / df

  result <- list(sigma.g2.hat = sigma.g2.hat, d0 = df2.corrected, s02 = s02)

  return(result)
}
