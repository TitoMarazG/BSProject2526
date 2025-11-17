P_DAG <- function(X, DAG, a, U) {
  
 #install.packages("bnlearn")
  source(Function_p(XA))
  
  j=nodes(DAG)
  
  i <- 1:length(j)
  pa <- parents(DAG,j[i])
  fa <- c(pa,j[i])  
  
  u <- 1:length(pa)
  v <- 1:length(fa)
  
  prod_pa <- prod(P_XA(X,a,U,pa[u]))
  prod_fa <- prod(P_XA(X,a,U,fa[v]))
  
  P_DAG = prod_fa/prod_p12a
  return(P_DAG)
                  
}
