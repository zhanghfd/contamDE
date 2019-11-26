##################################################### test between two conditions ########################################################
# tumor versus normal #
# tumor subtype1 versus tumor subtype2 #
contamDE.lm <-
function(counts, subtype=NULL, covariate=NULL, is.contaminated=TRUE, robust = TRUE) {
  counts <- as.matrix(counts)
  N <- ncol(counts) / 2  # the number of paired samples
  G <- nrow(counts)      # the number of total genes

  # use voom strategy in LIMMA package to estimate the initial value of purity proportion
  d <- limma_voom(counts = counts)
  p.limma <- d$p.limma
  log2FC.limma <- d$log2FC.limma

  ## design matrix ##
  # only one tumor subtype (without subtype information)
  if(is.null(subtype) | length(unique(subtype)) == 1) {
    subtype <- rep(1, N)
    n.sub <- 1
    if(is.null(covariate)) {
      design <- matrix(subtype, ncol = 1)
    } else {
      covariate <- matrix(covariate, ncol = N)
      design <- model.matrix(~ covariate)
    }
  } else {
  # with subtype information
    subtype <- factor(subtype)
    n.sub <- length(unique(subtype))
    if(is.null(covariate)){
      design <- model.matrix(~ 0 + subtype)
    } else {
      covariate0 <- matrix(covariate, nrow = N)
      design <- model.matrix(~ 0 + subtype + covariate0)
    }
  }

  K <- ncol(design)

  # sample size factor
  size <- d$size
  count.norm <- t(t(counts) / size)
  count.normal <- count.norm[,  (1: N)] + 1
  count.tumor  <- count.norm[, -(1: N)] + 1

  # do transformation
  y = log2(count.tumor) - log2(count.normal)

  w_hat <- rep(1,N)
  if(is.contaminated){
    ##########################################
    # proportion estimation
    ##########################################
    p.adj <- p.adjust(p.limma, method = 'fdr')
    log2FC <- log2FC.limma
    if(sum(p.adj < 0.1) > 1e3) p.adj[-order(p.adj)[1: 1e3]] <- 1

    up   <- which(p.adj < 0.1 & log2FC >  log2(1.5))
    down <- which(p.adj < 0.1 & log2FC < -log2(1.5))
    y_up <- y[up,]
    y_down <- y[down,]

    ## estimate ws
    sumup <- colSums(y_up)
    sumdown <- colSums(y_down)
    sum.max <- max(sumup - sumdown)
    w_hat <- (sumup - sumdown) / sum.max # purity proportion estimates
  }

  x_w <- design * w_hat
  dg <- N-K   # degree of freedom

  ## residual standard deviation
  res.ig <- (diag(N) - x_w %*% solve((t(x_w) %*% x_w)) %*% t(x_w)) %*% t(y)
  var.res <- colSums(res.ig ^ 2) / dg

  log_vg <- log(var.res) # log residual deviations

  mCoverage <- (count.normal + count.tumor) / 2;
  logC <- log(mCoverage)
  mlogC <- rowMeans(logC)

  ## gam fit
  gamres <- gam(y ~ s(x, k = 4), data = data.frame(x = mlogC, y = log_vg))
  predictors <- apply(logC, 2,function(x) {predict(gamres, data.frame(x = x))})

  v_ig <- exp(predictors)

  weight <- 1 / (v_ig) # the inverse of variances are the precision weights
  weight <- weight / rowMeans(weight)

  ## ebayes estimation

  ebayes.result <- ebayes(var.res = var.res, dg, robust = robust)
  hatsg2 <- ebayes.result$sigma.g2.hat
  df <- ebayes.result$d0 + dg

  ##t-test
  res <- matrix(0, G, (K + K ^ 2))
  for(g in 1: G) {
    y.g <- as.numeric(y[g, ])
    weight.gi <- weight[g, ]
    W.g <- diag(weight.gi)
    XWX.inv <- solve(t(x_w) %*% W.g %*% x_w)
    thetag.hat <- as.numeric(XWX.inv %*%t (x_w) %*% W.g %*% y.g)
    cov.theta <- XWX.inv * hatsg2[g]
    res[g, ] <- c(thetag.hat, cov.theta)
  }

  ##########################
  ##########################
  diag.id <- (0: (K - 1)) * K + (1: K)
  log2FC <- matrix(res[, (1: K)], nrow = G)
  log2FC.cov <- matrix(res[, -(1: K)], nrow = G)
  T.stat <- abs(log2FC) / sqrt(log2FC.cov[, diag.id])
  ##########################
  ##########################
  p.EB <- matrix(0, G, K)

  for(i in 1: K) {
    p.EB[, i] <- 2 * (1 - pt(T.stat[, i], df = df))
    large.T <- which(T.stat[, i] > 10)
    if(length(large.T) > 0)
      p.EB[large.T, i] <- 2 * (-pt(T.stat[large.T, i], df = df[large.T],log.p = TRUE))
  }

  if(is.null(covariate)) {
    colnames(p.EB) <- paste0('subtype-', 1:n.sub,' vs. normal')
  } else {
    colnames(p.EB) <- c(paste0('subtype-',   1: n.sub, ' vs. normal'),
                        paste0('covariate-', 1: ncol(covariate)))
  }
  d$p.contamDE.robust <- p.EB

############################################################
d$proportion <- w_hat
d$design <- design
d$log2FC <- log2FC
d$log2FC.cov <- log2FC.cov
d$df <- df
d$y <- y
d$weight <- weight
return(list(counts=d$counts,p.contamDE.lm=d$p.contamDE.robust,log2FC=d$log2FC,log2FC.cov=d$log2FC.cov,proportion=d$proportion,design=d$design,df=d$df,weight=d$weight,y=d$y))
}
