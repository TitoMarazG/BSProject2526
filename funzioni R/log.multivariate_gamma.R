log.multivariate_gamma <- function(p, x){
  val = 0
  for(j in 1:p)
    val = val + lgamma(x+(1-j)/2)
  return(p*(p-1)/4*log(pi)
         + val)
}