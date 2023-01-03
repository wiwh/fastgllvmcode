gradient_full <- function (fg) {
  # browser()
  # Compute Z
  fg$Z <- compute_Z(fg, start=fg$Z, maxit=10)$Z # TODO: maybe reduce to a maxit of 1?
  # cat(as.vector(cov(fg$Z)))
  # TODO: compare with the linpar obtained from comp_Z: if it is the same, take it.... and take the mean too
  # Update fg$linpar and fg$mean
  fg <- compute_mean(fg, linpar=NULL, return_object=T)

  dAB <- compute_dAB(fg, method="full")
  dphi <- compute_phi(fg)
  dcovZ <- cov(fg$Z)


  list(AB=dAB, phi=dphi, covZ = dcovZ)
}

gradient_simple <- function(fg) {
  # Update main values and compute gradient of the sample
  fg$Z <- compute_Z(fg, start=fg$Z, maxit=10)$Z # TODO: maybe reduce to a maxit of 1?
  # TODO: compare with the linpar obtained from comp_Z: if it is the same, take it.... and take the mean too
  fg <- compute_mean(fg, return_object=T) # TODO: check if we need this here
  # fg <- rescale(fg, rescale.A=T, rescale.B=T, target.cov=fg$parameters$covZ)

  dAB <- compute_dAB(fg, method="simple")
  dphi <- compute_phi(fg)
  dcovZ <- cov(fg$Z)

  list(AB=dAB, phi=dphi, covZ = dcovZ)
}


if(0) {
  devtools::load_all()
  set.seed(1234)
  poisson  <- 10
  gaussian <- 0
  binomial <- 0
  q <- 1
  k <- 0
  p <- poisson + gaussian + binomial
  family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial))
  set.seed(1030)
  fg <- gen_gllvmprime(nobs=100000, p=p, q=q, family=family, phi=3*(1:p)/p, k=k, A=matrix(sample(seq(-2, 2, l=p*q)),p, q), intercept=F, miss.prob = 0, scale=1)
  gradient <- gradient_full(fg)

  fg$parameters$A <- fg$parameters$A * 2
  set.seed(1030)
  fg <- simulate(fg, return_object = T)
  gradient2 <- gradient_full(fg)

  fg$parameters$A <- fg$parameters$A /4
  gradient3 <- gradient_full(fg)

  sim <- sapply(1:100, function(i){
    set.seed(i)
    if(i%%10==0 || i==1) cat("\n", i)
    fg_sim1 <- gen_gllvmprime(nobs=1000, p=p, q=q, A= fg$parameters$A, B = fg$parameters$B, phi=fg$parameters$phi, family=family, k=1, intercept=T)
    g1 <- gradient_full(fg_sim1)


    fg_sim2 <- fg_sim1
    fg_sim2$parameters$A <- fg_sim2$parameters$A *.9
    # fg_sim2 <- simulate(fg_sim2, return_oject = T)
    g2 <- gradient_full(fg_sim2)
    diff <- sapply(seq_along(g1), function(i) {
      g1[[i]] - g2[[i]]
    }, simplify=F)
    names(diff) <- names(g1)
    diff
  }, simplify=F)

  sim <- sapply(names(sim[[1]]),
                function(parname) do.call(rbind, sapply(sim, function(simi) as.vector(simi[[parname]]),
                                                        simplify=F)), simplify=F)

  par(mfrow=c(2,1))
  boxplot(sim$AB[,(order(fg$parameters$A))], outline=F); abline(h=0, col=2)
  boxplot(sim$phi);abline(h=0, col=2)
  par(mfrow=c(1,1))
}


