compute_dAB <- function(fg, method) {
  with(fg, {
    if (dimensions$k == 0) {
      ZX <- Z
    } else {
      ZX <- cbind(Z, X)
    }

    if ( method == "simple" ) {
      dAB <- t(Y) %*% (ZX/dimensions$n)
    } else if (method == "full") {
      dAB <- t(Y - mean) %*% (ZX/dimensions$n)
    } else {
      stop ("Unkown method `method`")
    }
    dAB
  })
}


if(0) {
  devtools::load_all()
  set.seed(1234)
  poisson  <- 0
  gaussian <- 10
  binomial <- 0
  q <- 2
  k <- 1
  p <- poisson + gaussian + binomial
  family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial))
  set.seed(1030)
  fg <- gen_gllvmprime(nobs=100, p=p, q=q, family=family, phi=3*(1:p)/p, k=k, intercept=T, miss.prob = 0, scale=1)
  psi <- compute_psi_simple(fg$Y, fg$X, fg$Z, fg$parameters, fg$families)
  gradient <- compute_gradients_simple(fg$Y, fg$X, fg$parameters, fg$families)

  # this must have expectaiton 0 under the true model. check this
  sim <- sapply(1:1000, function(i){
    set.seed(i)
    if(i%%10==0 || i==1) cat("\n", i)
    fg <- gen_gllvmprime(nobs=100, p=p, q=q, A= fg$parameters$A, B = fg$parameters$B, phi=fg$parameters$phi, family=family, k=1, intercept=T, miss.prob = 0, scale=1)
    psi <- compute_gradients_simple(fg$Y, fg$X, fg$parameters, fg$families, fg$Miss)
    psi
  }, simplify=F)

  sim <- sapply(names(sim[[1]]),
    function(parname) do.call(rbind, sapply(sim, function(simi) as.vector(simi[[parname]]),
    simplify=F)), simplify=F)

  image(cov(cbind(sim$A, sim$B))) # super big

  boxplot(sim$A, outline=F); abline(h=0, col=2)
  boxplot(sim$B,outline=F); abline(h=0, col=2)
  points(colMeans(sim$B), col=2)
  boxplot(sim$phi);abline(h=0, col=2)
}

if(0) {
  # Testing behavior of the different scores
  devtools::load_all()
  set.seed(1234)
  poisson  <- 100
  gaussian <- 0
  binomial <- 0
  q <- 1
  k <- 1
  p <- poisson + gaussian + binomial
  family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial))
  set.seed(12395)
  fg <- gen_gllvmprime(nobs=100, p=p, q=q, family=family, phi=3*(1:p)/p, k=k, intercept=T, miss.prob = 0, scale=1)

  fg$hessian <- simulate_hessian_AB(fg)
  for(i in 1:100){
    fg$hessian <- update_hessian_AB(fg$hessian, simulate_hessian_AB(fg), weight=.95)
  }
  # fg$hessian <- lapply(hessian, function(na) diag(1, ncol(na)))

  method <- "full"
  use_signs <- T
  hessian <- NULL

  set.seed(123)
  sims <- t(sapply(1:1000, function(na) {
    fg <- simulate(fg, return_gllvmprime = TRUE)
    controls <- list(method=method, use_signs=use_signs)
    as.vector(unlist(compute_dAB_centered(fg, controls=controls, hessian=hessian)[c("dA", "dB")]))
  }))
  {
  par(mfrow=c(2,1))
  boxplot(sims, outline=T, main=sqrt(mean(colMeans(sims)^2)))
  points(colMeans(sims), col=2)
  abline(h=0,col=2, lty=2)
  boxplot(sims, outline=F)
  points(colMeans(sims), col=2)
  abline(h=0,col=2, lty=2)
  par(mfrow=c(1,1))
  }
}
