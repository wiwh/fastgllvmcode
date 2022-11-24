devtools::load_all()
library(mirtjml)

rotate <- function(A, A.target){
  A.rot <- psych::Procrustes(A, A.target)$loadings
  error <- norm(A.rot - A.target)/norm(A.target)
  list(A=A.rot, error=error)
}

family <- binomial()

# from a big dataset of size n=50'000, p=2'500, we choose rho=(0.01, 0.02, ... 0.05) and compute the time to estimate
# a subset of the dataset of size ns=rho*n, ps=rho*p (rounded up).
# We compare est.sp, est.sa, est.gmf, est.gllvm

# We compare 2 methods:

est.mirtjml <- function(dat){
  tm <- proc.time()
  m <- mirtjml_expr(dat$Y, K=dat$dimensions$q)
  tmdiff <- proc.time()-tm

  rot <- rotate(m$A_hat, dat$parameters$A)

  list(A=rot$A, error=rot$error, time=tmdiff)
}

est.fastgllvm <- function(dat){
  tm <- proc.time()
  fit.sa <- fastgllvm(dat$Y, q = dat$dimensions$q, family="binomial",  intercept = T, hist=T, controls = list(maxit=100, alpha=1, beta=0, eps=1e-10, learning_rate.args=list(end=0.01, method="spall", rate=10)), method="simple", median=.2)
  tmdiff <- proc.time()-tm

  rot <- rotate(fit.sa$parameters$A, dat$parameters$A)

  list(A=rot$A, error=rot$error, time=tmdiff)
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
  list(fastgllvm=est.fastgllvm(dat),
       mirtjml = est.mirtjml(dat))
}

gen_dat <- function(n, p, rho){
  p.new <- ceiling(rho*p)
  n.new <- ceiling(rho*n)
  dat <- gen_fastgllvm(n.new,family="binomial", p=p.new, q=5, k=1, scale=.8)
  dat
}

onesim <- function(rho, n, p){
  library(mirtjml)
  devtools::load_all()
  set.seed(21414)
  dat <- gen_dat(n, p, rho)
  c(dat=dat, est.all(dat), nfull=n, pfull=p, rho, scale=0.8)
}

library(parallel)
# Detect the number of available cores and create cluster
cl <- parallel::makeCluster(detectCores()-2)
# Run parallel computation
parallel::clusterExport(cl, ls())

rhos <- seq(0.05, 1, by=0.05)
bigsim <- parallel::parLapply(cl,
                               rhos,
                               onesim, n=10000, p=2000)
save(bigsim, file="./experiments/mirtjlm/portion_data_results_large.Rdata")
parallel::stopCluster(cl)
