library(hdi)

library(pcalg)
library(Rgraphviz)
source("C:/Users/arian/Desktop/Magistrale/Primo anno/Bayesian Statistics/Progetto/personalized_my_skeleton_senza_stable_fast.R")

data(riboflavin)

X <- riboflavin$x   # matrice delle covariate
y <- riboflavin$y   # vettore delle risposte

Xsc <- scale(as.matrix(X))
ysc <- scale(as.numeric(y))

q = ncol(X)
alpha_bayes = 0.11 # BF associated to posterior probability of edge inclusion of 0.9 [given prior probability of edge inclusion of 0.5]

skltn.bayes <- my_skeleton(X = Xsc, a = q, U = diag(1, q), 
                           alpha = alpha_bayes, 
                           BayesTest = BF_Gaussian,
                           saveBF = FALSE, verbose = FALSE)
#plot(skltn.bayes@graph)