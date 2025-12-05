BF <- function(XX, n, a, U, u, v, S){
  source("log.P_DAGM_puntuale.R")
  return (exp(log.P_DAGM_puntuale(XX, n, a, U, v, S) - log.P_DAGM_puntuale(XX,n, a, U, v, c(S, u))))
  #BF01 (a numeratore dag senza edge)
    
}
 
