library(gmf)
devtools::load_all()

n <- 100
p <- 500
q <- 2
k=0
set.seed(123)
dat <- gen_gllvm(n,family="bernoulli", p=p, q=q, k=k, scale=.5)
hist(sigmoid(dat$natpar))
if(0){
  dat <- gmf.simulate(n=n, m=p, d=k, p=q, family=binomial())
  dat$par$A <- dat$V
  dat$k <- k
  dat$q <- q
  dat$p <- p
}

#using Stochastic Approximation

par(mfrow=c(3,1))
ptm <- proc.time()
fit.ffa <- bernoulli.estimate.ffa(dat$Y, q=q, batch=ifelse(500/n<1, 500/n, F), X=dat$X, iter=200, reps=4, reps.decreased_rate=0.8,learning_rate.start=20, learning_rate.end=1, tol=1e-5)
time.ffa <- proc.time() - ptm

(proc.ffa <-procrustes(fit.ffa$A, dat$par$A))
ts.plot(fit.ffa$A.hist[,seq(1, p*q, l=100)])
ts.plot(fit.ffa$B.hist)

plot(dat$par$A, proc.ffa$Lrot)
abline(0,1,col=2)

# Using Sample Path
fit.sp <- bernoulli.estimate.sp(dat$Y, q, dat$X, H=10, reps=2, iter=100, learning_rate.start=20, learning_rate.end=1)
(proc.sp <-procrustes(fit.sp$A, dat$par$A))

ts.plot(fit.sp$A.hist[,seq(1, p*q, l=100)])
ts.plot(fit.sp$B.hist)
par(mfrow=c(1,1))
plot(dat$par$A, proc.sp$Lrot); abline(0,1,col=2)



if(0){

par <- dat$par
par$A <- fit$A
par$B <- fit$B
pred <- bernoulli.predict(dat$Y, par, verbose=T, eps=-1)
plot(pred$linpar, dat$natpar); abline(0,1,col=2)
plot(pred$Z, dat$Z)
}
# GMF
family <- binomial()
model.newton = gmf(Y = dat$Y, X = dat$X, d = dat$k, p=dat$q, maxIter = 1000,
                   family=family, method = "quasi", tol = 1e-5, intercept = F)
M.gmf.newton = model.newton$fit
(gmf.proc <- procrustes(model.newton$v, dat$par$A))

# gllvm
fit.gllvm <- gllvm(y= dat$Y, X= dat$X, formula = ~., num.lv=dat$q, family=binomial(link="logit"))
# fit.gllvm <- gllvm(y= dat$Y, num.lv=dat$q, family=binomial(link="logit"))
coefs <- coefficients(fit.gllvm)
(gllvm.proc <- procrustes(coefs$theta, dat$par$A))


plot(dat$par$A, proc.sp$Lrot);abline(0,1)
points(dat$par$A, gmf.proc$Lrot, col=2)
points(dat$par$A, gllvm.proc$Lrot, col=3)
