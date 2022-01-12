devtools::load_all()
library(gmf)
library(gllvm)

family <- binomial()


# for each sim_num in bigsim, there are 10 simsetting n. For each of these, there are


# from a big dataset of size n=50'000, p=2'500, we choose rho=(0.01, 0.02, ... 0.05) and compute the time to estimate
# a subset of the dataset of size ns=rho*n, ps=rho*p (rounded up).
# We compare est.sp, est.sa, est.gmf, est.gllvm

# We compare 4 methods:

est.gmf <- function(dat){
  set.seed(12451)
  tm <- proc.time()
  model.newton = gmf(Y = dat$Y, X = dat$X, d = dat$k, p=dat$q, maxIter = 250,
                     family=family, method="newton",gamma=1e-1, tol=1e-5, intercept = F)
  tmdiff <- proc.time()-tm

  A.est <- model.newton$v

  (proc <- procrustes(A.est, dat$par$A))
  list(A=proc$Lrot, B=model.newton$beta, error=proc$Error, time=tmdiff, deviance = model.newton$deviance)
}

# est.gllvm <- function(dat){
#   set.seed(12451)
#   tm <- proc.time()
#   fit.gllvm <- gllvm(y= dat$Y, X= dat$X, formula = ~., num.lv=dat$q, family=binomial(link="logit"))
#   tmdiff <- proc.time()-tm
#
#   coefs <- coefficients(fit.gllvm)
#   A.est <- coefs$theta
#
#   (proc <- procrustes(A.est, dat$par$A))
#   list(A=proc$Lrot, error=proc$Error, time=tmdiff)
# }

est.sp <- function(dat){
  set.seed(12451)
  tm <- proc.time()
  fit.sp <- bernoulli.estimate.sp(dat$Y, dat$q, dat$X, H=10, iter=250, learning_rate.start=100, learning_rate.end=10, verbose=T)
  tmdiff <- proc.time()-tm

  A.est <- fit.sp$A

  # ts.plot(fit.sp$A.hist)
  (proc <- procrustes(A.est, dat$par$A))
  # plot(dat$par$A, proc$Lrot); abline(0,1,col=2)
  list(A=proc$Lrot, error=proc$Error, time=tmdiff)
}

est.sa <- function(dat){
  set.seed(12451)
  tm <- proc.time()
  # perform a first run to provide a warm start
  fit.sa <- bernoulli.estimate.ffa(dat$Y, q=dat$q, batch=.1, X=dat$X, iter=50, reps=1, reps.decreased_rate=0.5,learning_rate.start=40, learning_rate.end=20, tol=1e-5, verbose=T)
  A.init <- fit.sa$A
  fit.sa <- bernoulli.estimate.ffa(dat$Y, q=dat$q, batch=1, X=dat$X, iter=200, reps=1, A.init=A.init, reps.decreased_rate=0.5,learning_rate.start=40, learning_rate.end=.0001, tol=1e-5, verbose=T)
  tmdiff <- proc.time()-tm


  # polyak averaging
  # fit.sa <- bernoulli.estimate.ffa(dat$Y, q=dat$q, batch=1, X=dat$X, iter=200, reps=1, A.init=A.init, reps.decreased_rate=0.5,learning_rate.start=40, learning_rate.end=10, tol=1e-5, verbose=T)
  # A.est <- matrix(Means(fit.sa$A.hist[-(1:100),]), p, q)

  A.est <- fit.sa$A

  # ts.plot(fit.sa$A.hist)
  (proc <- procrustes(A.est, dat$par$A))
  # plot(dat$par$A, proc$Lrot); abline(0,1,col=2)
  list(A=proc$Lrot, B=fit.sa$B, error=proc$Error, time=tmdiff)
}


gen_dat <- function(n, p, q, k=0){
   # par$A <- matrix(runif(p*q, -2, 2), p, q)
  dat <- gen_gllvm(n=n,family="bernoulli", par=par, scale=.8)
  dat
}

est.all <- function(dat){
  list(sa=est.sa(dat)
       # sp=est.sp(dat),
       # gllvm=est.gllvm(dat),
       #gmf=est.gmf(dat)
       )
}

onesim <- function(rho){
  library(gmf)
  library(gllvm)
  devtools::load_all()
  dat <- gen_dat(n=n*rho, p=p, q=q, k=0)
  c(dat=dat, est.all(dat), n=n*rho, p=p, q=q, rho=rho, scale=1)
}

if(0){
  # This is just to test whether onesim is doing OK.
  sim <- onesim(.1)
  plot(sim$dat.par$A, sim$sa$A)
  points(sim$dat.par$A, sim$gmf$A, col=2)
}


# model parameters
n <- 10000
p <- 500
q <- 5
k <- 0
par <- gen_par(p, q, k=k, scale=0.8)
par$A <- matrix(sample(c(-2, -1,  0, 1, 2), p*q, replace=T, prob = c(.1,.2,.4,.3,.1)), p, q)


# setting parameters
times <- 20
rhos <- seq(0.1, 1, by=0.1)

library(parallel)
# Detect the number of available cores and create cluster
cl <- parallel::makeCluster(detectCores()-2)
# Run parallel computation
parallel::clusterExport(cl, ls())

bigsim <- parallel::parLapply(cl, 1:times, function(seed){
  set.seed(131+seed)
  lapply(rhos, onesim)
})
save(bigsim, file="./experiments/portion_data_fixed_loadings/portion_data_fixed_loadings_many_sa_1-20.Rdata")

times <- 20
rhos <- seq(0.1, 1, by=0.1)
bigsim <- parallel::parLapply(cl, 1:times, function(seed){
  set.seed(131+seed)
  lapply(rhos, onesim)
})
save(bigsim, file="./experiments/portion_data_fixed_loadings/portion_data_fixed_loadings_many_sa_21-40.Rdata")
times <- 20
rhos <- seq(0.1, 1, by=0.1)
bigsim <- parallel::parLapply(cl, 1:times, function(seed){
  set.seed(131+seed)
  lapply(rhos, onesim)
})
save(bigsim, file="./experiments/portion_data_fixed_loadings/portion_data_fixed_loadings_many_sa_41-60.Rdata")
times <- 20
rhos <- seq(0.1, 1, by=0.1)
bigsim <- parallel::parLapply(cl, 1:times, function(seed){
  set.seed(131+seed)
  lapply(rhos, onesim)
})
save(bigsim, file="./experiments/portion_data_fixed_loadings/portion_data_fixed_loadings_many_sa_61-80.Rdata")
times <- 20
rhos <- seq(0.1, 1, by=0.1)
bigsim <- parallel::parLapply(cl, 1:times, function(seed){
  set.seed(131+seed)
  lapply(rhos, onesim)
})
save(bigsim, file="./experiments/portion_data_fixed_loadings/portion_data_fixed_loadings_many_sa_81-100.Rdata")
parallel::stopCluster(cl)
