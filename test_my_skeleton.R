setwd("~/Documents/Bayesian/Project/BSProject2526")
library(pcalg)
library(BiocGenerics)
source("personalized_my_skeleton_senza_stable_fast.R")

q = 5# Fix number of nodes

# fix the adjacency matrix of a DAG; alternatively 
# use BCDAG::rDAG()
D0 = matrix(c(rep(0,q),
              rep(0, q),
              rep(c(1,1,0,0,0),3)), byrow = T, nrow = q, ncol = q)
D0 # <- rDAG()

## The true DAG is D1, i.e. the one having an edge u -> v

# Generate L,D, i.e. the parameters of the modified Cholesky decomposition of Sigma
L = D0*matrix(runif(q*q, 0.1, 1), q, q)*   #rimuovo giusto l'intorno di 0
  sample(c(-1, 1), size = q*q, replace = TRUE)
diag(L) = 1
L

D = diag(rep(1, q))

# Build Sigma = (L')^-1 * D * (L)^-1
Sigma = solve(t(L))%*%D%*%solve(L) 
Sigma
n = 200

#Simulate data!

X = mvtnorm::rmvnorm(n, rep(0, q), Sigma) 
X
a <- 5
U <- diag(1e, q)

pc()
pc.prof <- pc(suffStat = list(C = cor(X), n = nrow(X)), indepTest = gaussCItest,
                         p = ncol(X),
                         alpha = 0.01)
skeleton.test <- my_skeleton(X, q, U, alpha = 10)
par(mfrow = c(1,3))
plot(as(D0, "graphNEL"), main = "")
plot(pc.prof, main = "")
plot(skeleton.test@graph, main = "Il mio Grafo Bayesiano")
dev.off()

# to-do implementation:
# ANDARE AVANTI COL REPORT
# salvare i BF dei test fatti (list molto lunga attenzione)
# qualche test/simulation

# per generare un dag:
# in bcdag
rDAG <- function (q, w) 
{
  DAG = matrix(0, q, q)
  colnames(DAG) = rownames(DAG) = 1:q
  DAG[lower.tri(DAG)] = stats::rbinom(n = q * (q - 1)/2, size = 1, 
                                      prob = w)
  return(DAG)
}
# esempio da far vedere: come al variare di n come migliora la stima del grafo
# misura di distanza tra grafi: structural hamming distance = shd() =
# = quante mosse "mancano" per arrivare da un grafo all'altro


# con calma provare a cimentarsi a convertirlo in C++




