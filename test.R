library(numDeriv)

gen_z <- function(n, q){
  matrix(rnorm(n*q), n , q)
}

gen_y  <- function(z, L, family=gaussian()){
  if(family$family=="gaussian"){
    generator <- function(linpar) apply(linpar, 2, function(linpar_j) rnorm(length(linpar_j), linpar_j, 1))
  } else if (family$family=="binomial"){
    generator <- function(linpar) apply(linpar, 2, function(linpar_j) rbinom(length(linpar_j),1,family$linkinv(linpar_j)))
  } else if (family$family=="poisson"){
    generator <- function(linpar) apply(linpar, 2, function(linpar_j) rpois(length(linpar_j),family$linkinv(linpar_j)))
  }
  linpar <- z %*% t(L)
  generator(linpar)
}

expect_one_summand <- function(L, family){
  z <- gen_z(n, q)
  y <- gen_y(z, L, family=family)
  # as.vector(t(y) %*% t(solve(t(L) %*% L + diag(q)) %*% t(L) %*% t(y)))
  as.vector(t(y - family$linkinv(z %*% t(L))) %*% z / nrow(z))
}

expect_prime_one_summand <- function(L, family){
  z <- gen_z(n, q)
  y <- gen_y(z, L, family=family)
  # as.vector(t(y) %*% t(solve(t(L) %*% L + diag(q)) %*% t(L) %*% t(y)))
  linpar <- z %*% t(L)
  as.vector(t(y) %*% z / nrow(z))
}

expect_prime <- function(L, family){
  family$mu.eta()
}

get_z <- function(y, L, family=gaussian()){
 t(apply(y, 1, function(yi)glm(yi~0+L, family=family$family)$coef))
}
get_z_2 <- function(y, L){
 t(solve(t(L) %*% L + diag(q)) %*% t(L) %*% t(y))
}

expect <- function(L, H=100, family=gaussian()){
  rowMeans(sapply(1:H, function(h){
    expect_one_summand(L, family=family)
  }))
}

n <- 1000
p <- 100
q <- 10
family <- binomial()
family <- gaussian()
family <- poisson()

jac <- function(theta, family=gaussian()){
  set.seed(123123)
  L <- matrix(theta, p,q)
  expect(L, family=family)
}

compute_hessian <- function(L, zi){
  linpar <-  as.vector(zi %*% t(L))
  t(L) %*% diag(family$mu.eta(linpar)) %*% L + diag(q)/p
}

set.seed(1231232)
L <- matrix(rnorm(p*q), p, q)
z <- gen_z(n, q)
y <- gen_y(z, L, family=family)


if(0){
  set.seed(123)
  N <- expect(L*2, family=family)
  set.seed(123)
  D <- expect(L, family=family)
  N/D/as.vector(L)

  zhat <- get_z(y, L, family=family)
  zhat2 <- get_z_2(y, L)

  plot(z, zhat)
  plot(z, zhat2, col=2)
  abline(0,1,col=2)

  par(mfrow=c(2,1))
  qqnorm(zhat)
  abline(0,1,col=2)
  qqnorm(zhat2)
  abline(0,1,col=2)
  par(mfrow=c(1,1))

  J <- jacobian(jac, as.vector(L), method.args=list(eps=1e-2), family=poisson())
  J <- jacobian(jac, as.vector(L), family=poisson())
  vec <- c(0.01, rep(0, length(L)-1))
  (jac(as.vector(L+vec), family=poisson())-jac(as.vector(L-vec), family=poisson()))/(2*vec[1])


  # Estimating Z


  i <- 2
  zi <- z[i,]
  yi <- y[i,]
  H <- compute_hessian(L, zi)
  zi <- zi + as.vector(solve(H) %*% (zi + t(L) %*% (yi-family$linkinv(as.vector(zi %*% t(L))))))

  plot(z[i,], zi)
}

zhat <- z
linpar <- zhat %*% t(L)
hess <- t(sapply(1:n, function(i){
  diag(solve(compute_hessian(L, zhat[i,])))
}))

i <- 2
yi <- y[i,]

zi <- z[i,]
zi <- rep(0, q)
zi.old <- zi
H <- compute_hessian(L, zi)
zi <- zi - as.vector(solve(H) %*% (t(L) %*% (yi-family$linkinv(as.vector(zi %*% t(L))))))

plot(z[i,], zi, main=sum(abs(zi.old-zi)))

zhat <- zhat + hess * (-zhat/p + (y-family$linkinv(linpar)) %*% L)
plot(z, zhat)

zhat2 <- get_z(y, L, family=family)
zhat <- zhat2

qn <- qnorm(1:n/(n+1))
z2 <- sapply(1:q, function(k){
  qn[order(order((y %*% L)[,k]))]
})

qn <- qnorm(1:(n*q)/(n*q+1))
c <- qn[order(order(y%*%L))]

plot(z, c)
points(z, zhat2, col=2)
