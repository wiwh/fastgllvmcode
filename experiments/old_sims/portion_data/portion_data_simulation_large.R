devtools::load_all()
library(gmf)
library(gllvm)

family <- binomial()

# from a big dataset of size n=50'000, p=2'500, we choose rho=(0.01, 0.02, ... 0.05) and compute the time to estimate
# a subset of the dataset of size ns=rho*n, ps=rho*p (rounded up).
# We compare est.sp, est.sa, est.gmf, est.gllvm

# We compare 4 methods:

est.gmf <- function(dat){
  tm <- proc.time()
  model.newton = gmf(Y = dat$Y, X = dat$X, d = dat$k, p=dat$q, maxIter = 500,
                     family=family, method="newton",gamma=1e-1, tol=1e-5, intercept = F)
  tmdiff <- proc.time()-tm

  A.est <- model.newton$v

  proc <- procrustes(A.est, dat$par$A)
  list(A=proc$Lrot, error=proc$Error, time=tmdiff, deviance = model.newton$deviance)
}

est.gllvm <- function(dat){
  tm <- proc.time()
  fit.gllvm <- gllvm(y= dat$Y, X= dat$X, formula = ~., num.lv=dat$q, family=binomial(link="logit"))
  tmdiff <- proc.time()-tm

  coefs <- coefficients(fit.gllvm)
  A.est <- coefs$theta

  proc <- procrustes(A.est, dat$par$A)
  list(A=proc$Lrot, error=proc$Error, time=tmdiff)
}

est.sp <- function(dat){
  tm <- proc.time()
  fit.sp <- bernoulli.estimate.sp(dat$Y, dat$q, dat$X, H=1, reps=2, iter=200, learning_rate.start=40, learning_rate.end=.1, verbose=T)
  tmdiff <- proc.time()-tm

  A.est <- fit.sp$A

  proc <- procrustes(A.est, dat$par$A)
  # plot(dat$par$A, proc$Lrot); abline(0,1,col=2)
  list(A=proc$Lrot, error=proc$Error, time=tmdiff)
}

est.sa.batch <- function(dat){
  tm <- proc.time()
  fit.sa.batch <- bernoulli.estimate.ffa(dat$Y, q=dat$q, batch=.05, X=dat$X, iter=250, reps=2, reps.decreased_rate=0.5,learning_rate.start=40, learning_rate.end=.01, tol=1e-5, verbose=T)
  tmdiff <- proc.time()-tm

  A.est <- fit.sa.batch$A

  proc <- procrustes(A.est, dat$par$A)
  points(dat$par$A, proc$Lrot, col=3); abline(0,1,col=2)
  list(A=proc$Lrot, error=proc$Error, time=tmdiff)
}

est.sa <- function(dat){
  tm <- proc.time()
  fit.sa <- bernoulli.estimate.ffa(dat$Y, q=dat$q, batch=F, X=dat$X, iter=250, reps=2, reps.decreased_rate=0.5,learning_rate.start=40, learning_rate.end=1, tol=1e-5, verbose=T)
  tmdiff <- proc.time()-tm

  A.est <- fit.sa$A

  proc <- procrustes(A.est, dat$par$A)
  #plot(dat$par$A, proc$Lrot); abline(0,1,col=2)
  list(A=proc$Lrot, error=proc$Error, time=tmdiff)
}

# takes the first rho-proportion of observations in the dataset
get_subset_data <- function(dat, rho){
  p.new <- ceiling(rho*dat$p)
  n.new <- ceiling(rho*dat$n)
  dat$n <- dat$n.new
  dat$p <- dat$p.new
  dat$par <- NULL
  dat$Y <- dat$Y[1:n.new, 1:p.new]
  dat$X <- dat$X[1:n.new]
}

est.all <- function(dat){
  list(sa=est.sa(dat),
       sp=est.sp(dat),
       # gllvm=est.gllvm(dat),
       gmf=est.gmf(dat))
}

gen_dat <- function(n, p, q){
  dat <- gen_gllvm(n=n,family="bernoulli", p=p, q=q, k=1, scale=.8)
  dat
}

onesim <- function(rho){
  library(gmf)
  library(gllvm)
  devtools::load_all()
  set.seed(21414)
  n = 10000
  p = 2000
  q = 20
  dat <- gen_dat(n*rho, p*rho, max(2, q*rho))
  c(dat=dat, est.all(dat), n=n*rho, p=p*rho, q=min(2,q*rho), rho=rho, scale=0.8)
}

library(parallel)
# Detect the number of available cores and create cluster
cl <- parallel::makeCluster(detectCores()-2)
# Run parallel computation
parallel::clusterExport(cl, ls())

rhos <- seq(0.01, 0.095, by=0.005)
bigsim1 <- parallel::parLapply(cl,rhos, onesim)
save(bigsim1, file="./experiments/portion_data/portion_data_results_small_new2.Rdata")

# rhos <- seq(0.1, 1, by=0.05)
# bigsim2 <- parallel::parLapply(cl,
#                                    rhos,
#                                    onesim)
# save(bigsim2, file="./experiments/portion_data/portion_data_large_results_new.Rdata")
parallel::stopCluster(cl)
