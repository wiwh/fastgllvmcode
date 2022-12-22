devtools::load_all()
library(mirtjlm)
library(gmf)
library(gllvm)

family <- binomial()

# Simulation settings: (n, p, q, scale):
settings <- matrix(c(
  1000, 100, 5, .5,
  1000, 200, 5, .5,
  1000, 400, 5, .5,
  1000, 100, 10, .5,
  1000, 200, 10, .5,
  1000, 400, 10, .5
  ), ncol=4, byrow=T)

# We compare 4 methods:

est.gmf <- function(dat){
  tm <- proc.time()
  model.newton = gmf(Y = dat$Y, X = dat$X, d = dat$k, p=dat$q, maxIter = 500,
                     family=family, method="newton",  gamma=1e-2, tol = 1e-5, intercept = F)
  tmdiff <- proc.time()-tm

  A.est <- model.newton$v

  proc <- procrustes(A.est, dat$par$A)
  list(A=proc$Lrot, error=proc$Error, time=tmdiff)
}

est.sp <- function(dat){
  tm <- proc.time()
  fit.sp <- bernoulli.estimate.sp(dat$Y, q=dat$q, dat$X, H=10, reps=2, iter=100, learning_rate.start=40, learning_rate.end=1, verbose=T)
  tmdiff <- proc.time()-tm

  A.est <- fit.sp$A

  proc <- procrustes(A.est, dat$par$A)
  list(A=proc$Lrot, error=proc$Error, time=tmdiff)
}

est.sa <- function(dat){
  tm <- proc.time()
  fit.sa <- bernoulli.estimate.ffa(dat$Y, q=dat$q, batch=F, X=dat$X, iter=250, reps=2, reps.decreased_rate=0.5,learning_rate.start=40, learning_rate.end=1, learning_rate.type="exp", tol=1e-5, verbose=T)
  tmdiff <- proc.time()-tm

  A.est <- fit.sa$A

  proc <- procrustes(A.est, dat$par$A)
  plot(dat$par$A, proc$Lrot); abline(0,1,col=2)
  list(A=proc$Lrot, error=proc$Error, time=tmdiff)
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

est.sa.batch <- function(dat){
  tm <- proc.time()
  fit.sa <- bernoulli.estimate.ffa(dat$Y, q=dat$q, batch=ifelse(200<nrow(dat$Y), 200/nrow(dat$Y), F), X=dat$X, iter=250, reps=2, reps.decreased_rate=0.5,learning_rate.start=40, learning_rate.type="sa", learning_rate.end=.1, tol=1e-5, verbose=T)
  tmdiff <- proc.time()-tm

  A.est <- fit.sa$A

  proc <- procrustes(A.est, dat$par$A)
  list(A=proc$Lrot, error=proc$Error, time=tmdiff)
}

est.all <- function(dat){
  list(sa=est.sa(dat),
       #sp=est.sp(dat)),
       gmf=est.gmf(dat))
}

gen_dat <- function(n, p, q, scale){
  dat <- gen_gllvm(n,family="bernoulli", p=p, q=q, k=1, scale=scale)
  dat
}

onesim <- function(seed, n, p, q, scale=1){
  set.seed(seed)
  library(gmf)
  library(gllvm)
  devtools::load_all()
  dat <- gen_dat(n, p, q, scale=scale)
  c(dat=dat, est.all(dat), seed=seed, n=n, p=p, q=q, scale=scale)
}

library(parallel)
# Detect the number of available cores and create cluster
cl <- parallel::makeCluster(detectCores()-2)
# Run parallel computation
parallel::clusterExport(cl, ls())

bigsim <- vector(mode="list", l=nrow(settings))

tm <- proc.time()
for(i in seq_along(bigsim)){
  n <- as.numeric(settings[i,1])
  p <- as.numeric(settings[i,2])
  q <- as.numeric(settings[i,3])
  s <- as.numeric(settings[i,4])
  bigsim[[i]] <- parallel::parLapply(cl,
                      1:1,
                      onesim, n=n, p=p, q=q, scale=s)
  # Close cluster
}
tm <- proc.time()-tm
save(bigsim, file="./experiments/high_dimensions/high_dimensions_results.Rdata")
parallel::stopCluster(cl)
