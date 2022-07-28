devtools::load_all()
set.seed(123)
fg <- gen_fastgllvm(nobs=1000, p=300, q=1, family=c(rep("poisson", 300), rep("gaussian", 000), rep("binomial", 0)), k=1, intercept=F, miss.prob = 0)
zhat <- with(fg, compute_zstar(Y, parameters$A, parameters$phi, linpar$XB, families, Miss=Miss))
plot(fg$Z, zhat$Zstar)


set.seed(304434)
sim <- sapply(1:1000, function(aa){
  cat("\n", aa)
  Y_sim <- with(fg, {
    generate_y(
      linpar = NULL,
      phi = parameters$phi,
      families = families,
      A = parameters$A,
      B = parameters$B,
      X = X,
      Z = NULL,
      nobs = nrow(Y),
      Miss = Miss
    )
  })
  as.vector(compute_psi_star(Y_sim$Y, fg$X, Y_sim$Z, fg$parameters, fg$families, fg$Miss, compute_hessian=F)$psi_AB)
}, simplify=F)

sim <- do.call(rbind, sim)

plot(as.vector(fg$parameters$A), colMeans(sim[,1:300]))

plot(colMeans(sim)[order(fg$parameters$A)])
plot(colMeans(sim)[order(order(fg$parameters$A))])

hist(colMeans(sim), main=mean(sim), breaks=20)


#0.33, 0.30
#200: -0.01, -0.02
