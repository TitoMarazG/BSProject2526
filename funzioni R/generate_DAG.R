generate_DAG <- function(q, seed, w, n){
  # q: Numero nodi
  # w: Prob di un arco nel DAG generato casualmente 
  # n: numero dati simulati
  set.seed(seed)
  D0 <- rDAG(q, w)  # da BCDAG

  # Parametri della decomposizione di Cholesky modificata di Sigma
  # Matrice L (lower triangular con struttura del DAG)
  L <- D0 * # Matrice di adicenaza
    matrix(runif(q * q, 0.1, 1), q, q) *  # Coefficienti uniformi [0.1, 1]
    sample(c(-1, 1), size = q * q, replace = TRUE)  # Segni casuali
  diag(L) <- 1

  # Matrice D (diagonale)  
  D <- diag(rep(1, q)) 

  # Matrice di covarianza: Sigma = (L')^(-1) * D * L^(-1)
  Sigma <- chol2inv(t(L)) %*% D %*% chol2inv(L) 

  # Generazione dataset multivariato normale
  X <- mvtnorm::rmvnorm(n = n, mean = rep(0, q), sigma = Sigma)
  return(list(DAG_vero = D0,data = X))
}


  