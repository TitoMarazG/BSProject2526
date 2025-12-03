log.norm_cost <- function(a, U){
  source("log.multivariate_gamma.R")
  q = ncol(U)
  return(a/2*log(det(U))
         - a*q/2*log(2)
         - log.multivariate_gamma(q, a/2))
}