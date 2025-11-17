P_DAG <- function(X, DAG, a, U) {
  
 #install.packages("bnlearn")
  source(Function_p(XA))
  
  j=nodes(DAG)
  prod_pa = 1
  prod_fa = 1
  
 for(i in 1:length(j)){
   
   pa <- parents(DAG,j[i])
   fa <- c(pa,j[i])  
   
   prod_pa <-  prod_pa * (P_XA(X,a,U,pa))
   prod_fa <- prod_fa * (P_XA(X,a,U,fa))
   
 }
  
  P_DAG = prod_fa/prod_pa
  return(P_DAG)
                  
}
