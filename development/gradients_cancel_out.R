devtools::load_all()

set.seed(214)
p <- 100
q <- 5

par(mfrow=c(3,1))
for(n in c(100, 1000, 10000)){
  gllvm <- gen_gllvm(n=n, p=p, q=q, family="bernoulli")
  Y <- gllvm$Y
  X <- gllvm$X
  par <- gllvm$par

  grad <- gradient_bernoulli_A(Y, X, par)
  Egrad <- expected_gradient_bernoulli_A(X, par)

  plot(grad, Egrad, main=paste("n = ", n)); abline(a=0,b=1, col=2)
}
par(mfrow=c(1,1))
