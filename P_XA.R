p_XA <- function(X, a, U, A){
  n = nrow(X)
  q = ncol(X)
  A <- sort(as.integer(A)) # per evitare errori in chiamate come X[, A] o S[A, A]
  card_A = length(A)
  card_nonA = q - card_A
  S = t(X) %*% X # S = X'X
  
  S_AA = S[A, A, drop = FALSE]
  U_AA = U[A, A, drop = FALSE]
  
  c_prior <- norm_cost(a - card_nonA, U_AA)
  c_post  <- norm_cost(a + n - card_nonA, U_AA + S_AA)
  
  
  factor = (2 * pi)^(- (n * card_A) / 2)
  p_val = factor * (c_prior / c_post)
  return(p_val) 
}