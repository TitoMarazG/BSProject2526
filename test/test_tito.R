library(BCDAG)
library(pcalg)

source("my_skeleton.R")

# Simulate data with a given graphical structure --------------------------

# Adjacency matrix of a DAG (manual)

q = 5 # Number of nodes

D0 = matrix(c(rep(0,q),
              rep(0, q),
              rep(c(1,1,0,0,0),3)), byrow = T, nrow = q, ncol = q) 

# Adjacency matrix of a DAG (random)

prob = 0.3 # Probability of edge inclusion
D1 = rDAG(q, prob)

# Generate L,D, i.e. the parameters of the modified Cholesky decomposition of Sigma
adjM = D1

L = adjM*matrix(runif(q*q, 1, 3), q, q)*
  sample(c(-1, 1), size = q*q, replace = TRUE)
diag(L) = 1

D = diag(rep(1, q))

# Build Sigma = (L')^-1 * D * (L)^-1
Sigma = solve(t(L))%*%D%*%solve(L) 


# Data simulation
n = 200 # Sample size

X = mvtnorm::rmvnorm(n, rep(0, q), Sigma)


# PC Algorithm (frequentist) ----------------------------------------------
suffStat = list(C = cor(X), n = nrow(X))
indepTest = gaussCItest
alpha = 0.05

skeleton.freq = skeleton(suffStat = suffStat,
                             indepTest = indepTest,
                             p = ncol(X),
                             alpha = alpha, method = "stable.fast")

pc.freq = pc(suffStat = suffStat,
             indepTest = indepTest,
             p = ncol(X),
             alpha = alpha)

# Plots
par(mfrow = c(2,2))
plot(skeleton.freq, main = "Frequentist PC-estimated Skeleton")
plot(as(adjM, "graphNEL"), main = "True CPDAG")


plot(pc.freq, main = "Frequentist PC-estimated CPDAG")
plot(as(adjM, "graphNEL"), main = "True CPDAG")


# PC Algorithm (Bayesian) -------------------------------------------------
a = q
U = diag(1, q)
BF0 = 10

skeleton.bayes = my_skeleton(X, a, U, alpha = BF0, method = "original")






