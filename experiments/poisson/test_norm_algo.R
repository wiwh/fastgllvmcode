get_K <- function(A, Psi){
  p <- nrow(A)
  q <- ncol(A)
  solve(t(A) %*% diag(Psi**-1) %*% A + diag(q)) %*% t(A) %*% diag(Psi**-1)
}

get_Sigma <- function(A, Psi){
  A %*% t(A) + diag(Psi)
}

gen_A <- function(p, q){
  matrix(rnorm(p*q),p,q)
}

get_S <- function(p, heywood=F){
  n <- 5*p
  S <- matrix(rnorm(n*p), n, p)
  S <- t(S) %*% S /n
  if(heywood){
    es <- eigen(S)
    es$values[round(p*.9):p] <- 1e-4
    S <- es$vectors %*% diag(es$values) %*% t(es$vectors)
  }
  S
}

get_Psi <- function(S, A){
  Psi <- diag(S - A %*% t(A))
  Psi[Psi<1e-2] <- 1e-2
  Psi
}

Msqrt <- function(A){
  e <- eigen(A)
  if(any(e$values<1e-4)){
    warning("Eigenvalues smaller than 1e-4: smoothing was done.")
    e$values[e$values < 1e-2] <- 1e-2
  }
  Msqrt <- e$vectors %*% diag(e$values^(1/2)) %*% t(e$vectors)
  Msqrtneg <- e$vectors %*% diag(e$values^(-1/2)) %*% t(e$vectors)
  list(Msqrt = Msqrt, Msqrtneg=Msqrtneg)
}

LHS <- function(S, A1, A2, method="chol"){
  Psi1 <- get_Psi(S, A1)
  Psi2 <- get_Psi(S, A2)
  Sig1 <- get_Sigma(A1, Psi1)
  Sig2 <- get_Sigma(A2, Psi2)
  K1 <- get_K(A1, Psi1)
  K2 <- get_K(A2, Psi2)
  if(method=="chol"){
    A1next <- S %*% t(K1) %*% t(solve(chol(K1 %*% S %*% t(K1)))) %*% chol(K1 %*% Sig1 %*% t(K1))
    A2next <- S %*% t(K2) %*% t(solve(chol(K2 %*% S %*% t(K2)))) %*% chol(K2 %*% Sig2 %*% t(K2))
  } else if (method=="eigen"){
    S1.sqrt <- Msqrt(K1 %*% S %*% t(K1))
    S2.sqrt <- Msqrt(K2 %*% S %*% t(K2))
    Sig1.sqrt <- Msqrt(K1 %*% Sig1 %*% t(K1))
    Sig2.sqrt <- Msqrt(K2 %*% Sig2 %*% t(K2))
    A1next <- S %*% t(K1) %*% S1.sqrt$Msqrtneg %*% Sig1.sqrt$Msqrt
    A2next <- S %*% t(K2) %*% S2.sqrt$Msqrtneg %*% Sig2.sqrt$Msqrt
  }
  
  norm((A1next - A1) - (A2next -A2), type="F")
}

ffa.eigen <- function(S, q, method="eigen"){
  p <- nrow(S)
  A <- gen_A(p, q)
  # A <- diag(1, p, q)
  iter <- 500
  conv  <- rep(1, iter)
  for(i in 1:iter){
    A.old <- A
    Psi <- get_Psi(S, A)
    Sig <- get_Sigma(A, Psi)
    K <- get_K(A, Psi)
    if(method=="eigen"){
      S.sqrt <- Msqrt(K %*% S %*% t(K))
      Sig.sqrt <- Msqrt(K %*% Sig %*% t(K))
      A <- S %*% t(K) %*% S.sqrt$Msqrtneg %*% Sig.sqrt$Msqrt
    }
    if(method=="chol"){
      S.c <- chol(K %*% S %*% t(K))
      Sig.c <- chol(K %*% Sig %*% t(K))
      A <- S %*% t(K) %*% t(solve(S.c)) %*% Sig.c
    }
    conv[i] <- norm(A - A.old, type="F")
  }
  list(A=A, Psi=Psi, covmat = A %*% t(A) + diag(Psi))
}


psych.fa <- function(S, q, fm="minres"){
  fit <- psych::fa(S, nfactors = q, covar=T, fm=fm)
  list(A=fit$loadings, Psi = fit$uniquenesses, covmat = fit$loadings %*% t(fit$loadings) + diag(fit$uniquenesses))
}


ffa <- function(S, q){
   fit <- fastfactoranalysis::ffa(x = 1, nfactors = q, covx = S)
   A <- fit$loadings
   Psi <- fit$communalites
   list(A=A, Psi=Psi, covmat=A %*% t(A) + diag(Psi))
}


# set.seed(125)
p <- 20 
q <- 12 
S <- get_S(p, heywood=F)

# fac <- factanal(factors=q, covmat=S)
ffe <- ffa.eigen(S, q, method="eigen")
ff <- ffa(S, q)
ps <- psych.fa(S, q, fm="minres")

plot(S, ffe$covmat)
points(S, ps$covmat, col=2)
points(S, ff$covmat, col=3)


norm(S - ffe$covmat, type="F")
norm(S - ps$covmat, type="F")
norm(S - ff$covmat, type="F")

RHS <- function(S, A1, A2){
  Psi1 <- get_Psi(S, A1)
  Psi2 <- get_Psi(S, A2)
  K1 <- get_K(A1, Psi1)
  K2 <- get_K(A2, Psi2)
  norm(A1 - A2, type="F")
}


q <- 5 
p <- 20 

S <- get_S(p)
A <- gen_A(p, q)
Psi <- get_Psi(S, A)
K <- get_K(A, Psi)
Sig <- get_Sigma(A, Psi)


sims <- sapply(1:1000, function(i){
  set.seed(i)
  A1 <- gen_A(p, q) * runif(1,1,2)*100 # /1000 * rexp(1, 100)
  A2 <- gen_A(p, q) * runif(1,1,2)*100 #*1000 * rexp(1, 1/100)
  
  
  S <- get_S(p) * runif(1, 1,2)*1000 #* rexp(1, 1/100) #  + (diag(A1 %*% t(A1)) + diag(A2 %*% t(A2)))/5
  
  Psi1 <- get_Psi(S, A1)
  Psi2 <- get_Psi(S, A2)
  
  L1 <- LHS(S, A1, A2, method="eigen")
  R1 <- RHS(S, A1, A2)
  
  
  cat("\n", 123+i,  L1 < R1, round(L1/R1, 4))
  
  L1/R1
})

hist(sims, main=sum(sims>1))
