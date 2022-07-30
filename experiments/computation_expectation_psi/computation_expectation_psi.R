devtools::load_all()
set.seed(1423)
poisson  <- 20
gaussian <- 000
binomial <- 0
p <- poisson + gaussian + binomial
fg <- gen_fastgllvm(nobs=1000, p=p, q=2, family=c(rep("poisson", poisson), rep("gaussian", gaussian), rep("binomial", binomial)), k=1, intercept=F, miss.prob = 0)
# zhat <- with(fg, compute_zstar(Y, parameters$A, parameters$phi, linpar$XB, families, Miss=Miss))
# plot(fg$Z, zhat$Zstar)


factor <- 0.8
# we have a problem right here
par.init <- fg$parameters
par.init$A <- fg$parameters$A*factor

Zstar <- compute_zstar(fg$Y, par.init$A, fg$parameters$phi, fg$X %*% t(par.init$B), fg$families)$Zstar
plot(fg$Z, Zstar, ylim=c(-3,3), xlim=c(-3,3), main=factor); abline(0,1,col=2)

linpar_bprime <- compute_linpar_bprime(compute_linpar(Zstar, par.init$A, fg$X, par.init$B)$linpar, fg$families)

i <- i+1
set.seed(31231+i)
Y_sim <- generate_y(
    linpar = NULL,
    phi = par.init$phi,
    families = fg$families,
    A = par.init$A,
    B = par.init$B,
    X = fg$X,
    Z = NULL,
    nobs = nrow(fg$Y),
    Miss = fg$Miss
)
Zsim <- compute_zstar(Y_sim$Y, par.init$A, fg$parameters$phi, fg$X %*% t(par.init$B), fg$families, start=Y_sim$Z)$Zstar
linpar_bprime_sim <- compute_linpar_bprime(compute_linpar(Zsim, par.init$A, fg$X, par.init$B)$linpar, fg$families)

psi_est_hess <- compute_psi_star(fg$Y, fg$X, fg$Z, par.init, fg$families, fg$Miss, compute_hessian=T)

psi_tru <- t(fg$Y-fg$linpar$linpar) %*% cbind(fg$Z, fg$X)
psi_est <- t(fg$Y-linpar_bprime) %*% cbind(Zstar, fg$X)
psi_sim <- t(Y_sim$Y-linpar_bprime_sim) %*% cbind(Y_sim$Z, fg$X)

psi_est_hess <- compute_psi_star(Ysim_$Y, fg$X, Z=NULL, par.init, fg$families, fg$Miss, compute_hessian=T)
psi_est_hess <- compute_psi_star((fg$Y + Y_sim$Y)/2, fg$X, Z=NULL, par.init, fg$families, fg$Miss, compute_hessian=T)
hess_sam <- compute_psi_star(fg$Y, fg$X, Z=NULL, par.init, fg$families, fg$Miss, compute_hessian=T)
hess_est <- compute_psi_star(Y_sim$Y, fg$X, Z=NULL, par.init, fg$families, fg$Miss, compute_hessian=T)

lapply(seq_along(hess_sam$psi_AB_hessian))

psi_est <- t(fg$Y) %*% cbind(Zstar, fg$X)
psi_sim <- t(Y_sim$Y) %*% cbind(Y_sim$Z, fg$X)
plot(cbind(fg$parameters$A, fg$parameters$B), psi_est);abline(0,1,col=2)
plot(cbind(fg$parameters$A, fg$parameters$B), -compute_hessian_x_psi(psi_est - psi_sim, psi_est_hess$psi_AB_hessian), main=factor); abline(h=0,col=2)

# This works!

# Another idea: 1) estimate Z, 2) rescale Z and A 3) update A.

Zstar <- compute_zstar(fg$Y, par.init$A, fg$parameters$phi, fg$X %*% t(par.init$B), fg$families)$Zstar
plot(fg$Z, Zstar); abline(0,1,col=2)
res <- rescale(Zstar, par.init$A)
par.init$A <- res$A
Zstar <- res$Z
points(fg$Z, Zstar, col=2)
psi_est_hess <- compute_psi_star(fg$Y, fg$X, Z=Zstar, par.init, fg$families, fg$Miss, compute_hessian=T)
psi_weighted <- compute_hessian_x_psi(psi_est_hess$psi_AB, psi_est_hess$psi_AB_hessian)
plot(cbind(fg$parameters$A, fg$parameters$B), psi_est_hess$psi_AB)
plot(cbind(fg$parameters$A, fg$parameters$B), psi_weighted)


psi_tru <- t(fg$Y) %*% fg$Z
psi_est <- t(fg$Y) %*% Zstar

fg_parameters$


plot(linpar_bprime$tru, linpar_bprime$est, col=1, main=factor); abline(0,1,col=2)
plot(fg$linpar$linpar, linpar_bprime$est, col=1, main=factor); abline(0,1,col=2)

psi_sam <- compute_psi_star(fg$Y, fg$X, fg$Z, par.init, fg$families, fg$Miss, compute_hessian=T)
psi_sam_hess <- compute_hessian_x_psi(psi_sam$psi_AB, psi_sam$psi_AB_hessian)

par(mfrow=c(1,2))
plot(fg$parameters$A, psi_sam_hess[,1])
plot(fg$parameters$A, psi_sam$psi_AB[,1])
par(mfrow=c(1,1))


set.seed(304434)
sim <- sapply(1:1000, function(aa){
  cat("\n", aa)
  Y_sim <- with(fg, {
    generate_y(
      linpar = NULL,
      phi = par.init$phi,
      families = families,
      A = par.init$A,
      B = par.init$B,
      X = X,
      Z = NULL,
      nobs = nrow(Y),
      Miss = Miss
    )
  })
  # zhat <- compute_zstar(Y_sim$Y, par.init$A, par.init$phi, Y_sim$linpar$XB, fg$families, start=Y_sim$Z, Miss=fg$Miss)
  # plot(Y_sim$Z, zhat$Zstar); abline(0,1,col=2)
  psi_star <- compute_psi_star(Y_sim$Y, fg$X, Y_sim$Z, par.init, fg$families, fg$Miss, compute_hessian=T)
  as.vector(compute_hessian_x_psi(psi_star$psi_AB, psi_star$psi_AB_hessian))
}, simplify=F)
sim <- do.call(rbind, sim)

boxplot(sim[,order(par.init$A)]);abline(h=0,col=2)

plot(par.init$A, colMeans(sim[,1:p]))

sim_1.5_maxit_1 <- sim
plot(fg$parameters$A, colMeans(sim_1.5_maxit_1[,1:p]))
plot(fg$parameters$A, colMeans(sim_1.5[,1:p]))
plot(fg$parameters$A, colMeans(sim_1[,1:p]))
plot(fg$parameters$A, colMeans(sim_0.5[,1:p]))

plot(fg$parameters$B, colMeans(sim_1.5[,-(1:p)]))
plot(fg$parameters$B, colMeans(sim_1[,-(1:p)]))
plot(fg$parameters$B, colMeans(sim_0.5[,-(1:p)]))
