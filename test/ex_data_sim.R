# The goal of this script is to simulate multivariate normal data (q=5),
# which will be contained in the X matrix at the end (n=200 simulations),
# based on the DAG D0 defined at the beginning

library(BCDAG)
q = 5# Fix number of nodes

# fix the adjacency matrix of a DAG; alternatively 
# use BCDAG::rDAG()
D0 = matrix(c(rep(0,q),
              rep(0, q),
              rep(c(1,1,0,0,0),3)), byrow = T, nrow = q, ncol = q)
D0

## The true DAG is D1, i.e. the one having an edge u -> v

# Generate L,D, i.e. the parameters of the modified Cholesky decomposition of Sigma
L = D0*matrix(runif(q*q, 1, 3), q, q)*
  sample(c(-1, 1), size = q*q, replace = TRUE)
diag(L) = 1
L

D = diag(rep(1, q))

# Build Sigma = (L')^-1 * D * (L)^-1
Sigma = solve(t(L))%*%D%*%solve(L) 

n = 200

#Simulate data!

X = mvtnorm::rmvnorm(n, rep(0, q), Sigma) 
X
