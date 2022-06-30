gen_y  <- function(x, beta, family=gaussian()){
  if(family$family=="gaussian"){
    generator <- function(linpar) apply(linpar, 2, function(linpar_j) rnorm(length(linpar_j), linpar_j, 1))
  } else if (family$family=="binomial"){
    generator <- function(linpar) apply(linpar, 2, function(linpar_j) rbinom(length(linpar_j),1,family$linkinv(linpar_j)))
  } else if (family$family=="poisson"){
    generator <- function(linpar) apply(linpar, 2, function(linpar_j) rpois(length(linpar_j),family$linkinv(linpar_j)))
  }
  linpar <- x %*% beta
  generator(linpar)
}


# score and hessian function

score <- function(y, x, beta, family){
  t(x) %*% (y - family$linkinv(x %*% beta))  - beta
}
hessian <- function(y, x, beta, family){
  -t(x) %*% (x*as.vector(family$mu.eta(x %*% beta))) - diag(p)
}

newton <- function(y, x, start=start_values(y,x), family=gaussian(), maxit=100, diag=F, line.search=T, thresh=1e-3){
  beta <- start
  score.beta <- function(beta) score(y, x, beta, family)
  hist <- list()
  for(i in 1:maxit){
    # print(beta)
    beta.old <- beta
    myscore <- score.beta(beta)
    if(!diag | length(beta)==1){
      (beta <- beta - as.vector(solve(hessian(y, x, beta, family)) %*% myscore))
    }else{
      (beta <- beta - as.vector(diag(1/diag(hessian(y, x, beta, family))) %*% myscore))
    }

    # line search
    if(line.search) beta <- line_search(beta.old, beta, score.beta)

    hist[[i]] <- beta
    if(mean(abs(beta.old-beta))<thresh)break()
  }
  hist <- do.call(rbind, hist)
  list(beta=beta, hist=hist, niter=i)
}

start_values <- function(y, x, type=1){
  if(type==1)start <- as.vector(lm(log(y+1e-3)~0+x)$coef)
  start
}

line_search <- function(b1, b2, f){
  t <- c(0.5, 0.75, 1, 1.25, 1.5, 2, 5, 10)
  search <- sapply(t, function(ti){
    mean(f(b1 + ti*(b2-b1))^2)
  })
  cat("\n", t[which.min(search)])
  search[is.na(search)] <- Inf
  b1 + t[which.min(search)]*(b2-b1)
}

n <- 100
p <- 10
set.seed(1231)
x <- matrix(rnorm(n*p), n, p)*4
set.seed(32192)
beta0 <- rnorm(p)/10

family <- binomial()
family <- poisson()

set.seed(123)
y <- as.vector(gen_y(x, beta0, family))

(f1 <- newton(y,x,family=family, beta0, line.search = F, diag=F, thresh=1e-3))

library(glmnet)
(b1 <- as.vector(glmnet(x, y, family=family, alpha=0, lambda=1/n, thresh=1e-10)))
(b2 <- glm(y~0 +x, family=family)$coef) # , start=start_values(y,x,type=1)))
# (b3 <- newton(y,x, family=family, line.search = T))
ts.plot(f1$hist)

# plot(beta0, f1$beta)
# points(beta0, b2, col=2)
# abline(0,1,col=1)

stop()
(f2 <- newton(y,x,family=family, line.search = F, diag=F, maxit = 1))
(f3 <- newton(y,x,family=family,beta0, line.search = F, diag=T))


# (newton(y,x,beta0, family, line.search = F, diag=F, maxit=11990))

microbenchmark::microbenchmark(newton(y,x,beta0,family, line.search=T))
microbenchmark::microbenchmark(newton(y,x,beta0,family, line.search=F))


path <- sapply(seq(1, 1000, l=20), function(s){
  newton(y,x, family=family, start=beta0, maxit=s, diag=T)
})
ts.plot(t(path))

f1 <- function(x,y,family)as.vector(glmnet(x, y, family=family, alpha=0, lambda=1/n, thresh=1e-10)$beta)
f1.2 <- function(x,y,family)as.vector(glmnet(x, y, family=family, alpha=1, lambda=1/n, theresh=1e-10)$beta)
f2 <- function(x,y,family)as.vector(glmnet(x, y, family=family, alpha=0)$beta)
f3 <- function(x,y,family)as.vector(glm(y~0+x, family=family)$coef)
microbenchmark::microbenchmark(f1.2(x,y,family), times=10)

microbenchmark::microbenchmark(f1(x,y,family),f2(x,y,family), f3(x,y,family), times=10)

f1(x,y,family)
