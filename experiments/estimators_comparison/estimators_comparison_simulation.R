devtools::load_all()
library(gmf)
library(gllvm)

family <- binomial()

# Simulation settings:
rep <- 1000
q <- 2
p.list <- c(20, 40)
n.list <- c(50, 100, 500)
# We compare 4 methods:

est.gmf <- function(dat){
  tm <- proc.time()
  model.newton = gmf(Y = dat$Y, X = dat$X, d = dat$k, p=dat$q, maxIter = 1000,
                     family=family, method = "quasi", tol = 1e-5, intercept = F)
  tmdiff <- proc.time()-tm

  A.est <- model.newton$v

  proc <- procrustes(A.est, dat$par$A)
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

est.sa <- function(dat){
  tm <- proc.time()
  fit.sa <- bernoulli.estimate.ffa(dat$Y, q=dat$q, batch=F, X=dat$X, iter=200, reps=4, reps.decreased_rate=0.8,learning_rate.start=50, learning_rate.end=1, tol=1e-5, verbose=F)
  tmdiff <- proc.time()-tm

  A.est <- fit.sa$A

  proc <- procrustes(A.est, dat$par$A)
  list(A=proc$Lrot, error=proc$Error, time=tmdiff)
}

est.sp <- function(dat){
  tm <- proc.time()
  fit.sp <- bernoulli.estimate.sp(dat$Y, dat$q, dat$X, H=10, reps=2, iter=100, learning_rate.start=20, learning_rate.end=1, verbose=F)
  tmdiff <- proc.time()-tm

  A.est <- fit.sp$A

  proc <- procrustes(A.est, dat$par$A)
  list(A=proc$Lrot, error=proc$Error, time=tmdiff)
}

est.all <- function(dat){
  list(sp=est.sp(dat),
       sa=est.sa(dat),
       gmf=est.gmf(dat),
       gllvm=est.gllvm(dat))
}

gen_dat <- function(n, p){
  par <- gen_par(p,q, k=1)
  par$A[,1] <- seq(2,-2, l=p)
  par$A[,2] <- seq(-1,1, l=p)
  dat <- gen_gllvm(n,family="bernoulli", p=p, q=q, par=par)
  dat
}

onesim <- function(seed, n, p){
  set.seed(seed)
  library(gmf)
  library(gllvm)
  devtools::load_all()
  dat <- gen_dat(n, p)
  c(dat=dat, est.all(dat), seed=seed, n=n, p=p)
}

settings <- expand.grid(n.list, p.list)
colnames(settings) <- c("n", "p")

library(parallel)
# Detect the number of available cores and create cluster
cl <- parallel::makeCluster(detectCores()-2)
# Run parallel computation
parallel::clusterExport(cl, ls())

bigsim <- vector(mode="list", l=nrow(settings))

for(i in seq_along(bigsim)){
  n <- as.numeric(settings[i,1])
  p <- as.numeric(settings[i,2])
  bigsim[[i]] <- parallel::parLapply(cl,
                      1:500,
                      onesim, n=n, p=p)
  # Close cluster
}
save(bigsim, file="./experiments/estimators_comparison/estimators_comparison_results_setting2.Rdata")
parallel::stopCluster(cl)
