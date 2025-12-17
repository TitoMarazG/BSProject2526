BF <- function(XX, n, a, U, u, v, S){
  #XX = X(T) * X
  #n sample size
  #a prior degrees of freedom
  #U prior scale matrix
  #u nodo genitore
  #v nodo figlio (da testare)
  #S: parents di entrambi
  #D1: dag con edge in più u->v: viene formalizzato includendo u nel set S  
  source("log.P_DAGM_puntuale.R")
  return (exp(log.P_DAGM_puntuale(XX, n, a, U, v, S) - log.P_DAGM_puntuale(XX, n, a, U, v, c(S, u))))
  #BF01 (a numeratore dag senza edge)
}
 
