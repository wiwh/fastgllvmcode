library(gmf)
devtools::load_all()

n <- 1000
p <- 20
q <- 1
k=3
set.seed(123)
dat <- gen_gllvm(n,family="bernoulli", p=p, q=q, k=3, scale=1)

fit <- bernoulli.estimate.ffa(dat$Y, q, dat$X, iter=250)
fit2 <- bernoulli.estimate.sp(dat$Y, q, dat$X, H=2, reps=1, iter=250)

proc <-procrustes(fit$A, dat$par$A)
proc2 <-procrustes(fit2$A, dat$par$A)

plot(proc$Lrot, dat$par$A, main=proc$Error);abline(a=0,b=1,col=2);abline(0,-1,col=2)
plot(proc2$Lrot, dat$par$A, main=proc2$Error);abline(a=0,b=1,col=2);abline(0,-1,col=2)

par <- dat$par
par$A <- fit$A
par$B <- fit$B
pred <- bernoulli.predict(Y, par, verbose=T, eps=1e-15)

par <- dat$par
pred <- bernoulli.predict(Y, par, verbose=T, eps=1e-15)
plot(pred$linpar, dat$natpar)

# GMF
family <- binomial()
model.newton = gmf(Y = Y, X = dat$X, d = dat$k, p=q, gamma= 5e-1, maxIter = 1000,
                   family=family, method = "quasi", tol = 1e-5)
M.gmf.newton = model.newton$fit



# get goodness of fit
dev.null <- matrix.deviance(mean(Y), Y, family)
1 - matrix.deviance(sigmoid(pred$linpar), Y, family) / dev.null
1 - matrix.deviance(M.gmf.newton, Y, family) / dev.null




# In this example we compare a fitted model with a null model
data = gmf.simulate(family=binomial())

# Fit a GLLVM model using Newton method
model = gmf(data$Y, p=q)
fit <- bernoulli.ffa.covariate(data$Y, q=q,X=data$X)

# Compute mean deviance of the null model and the fitted model
dev.null = matrix.deviance(mean(data$Y), data$Y, model$family)
dev.model = matrix.deviance(model$fit, data$Y, model$family)

cat("Deviance explained by the model: ", 1-dev.model/dev.null)

isthisthefit <- apply(model$u %*% t(model$v), 1:2, sigmoid)
