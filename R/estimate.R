compute_K <- function(par){
  Ap <- A/Psi
  solve(t(A)%*% Ap + diag(1, q)) %*% t(Ap)
}

compute_learningRate <- function(iter, start, end){
  start*(exp(log(end/start)/iter))**(1:iter) - end
}

compute_gradient <- function(){

}

profile_score_gllvm <- function(){

}

score_gllvm <- function(params){

}
