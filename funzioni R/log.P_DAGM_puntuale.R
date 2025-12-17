log.P_DAGM_puntuale <- function(XX, n, a, U, v, S) {
  source("log.p_XA.R")
  return(log.p_XA(XX,n,a,U,c(S, v)) - log.p_XA(XX,n,a,U,S))
}
#calcola P(X|D) tramite = P(X|familyD(v))/P(X|parentsD(v))
#la D sta a significare che i parents/family sono presi dal DAG D
#XX = X(T) * X
#n = sample size
#a = hyperparameter
#U = hyperparameter
#u = nodo genitore
#v = nodo figlio (nodo da testare)
#S:parents di entrambi