# ==============================================================================
# SETUP INIZIALE
# ==============================================================================
setwd("~/Documents/Bayesian/Project/BSProject2526")

library(pcalg)
library(BiocGenerics)
library(BCDAG)

set.seed(123)  # Per riproducibilità

rm(list = ls()) 
graphics.off() 

source("/Users/leomarcellopoli/Documents/Bayesian/Project/BSProject2526/funzioni R/personalized_my_skeleton_senza_stable_fast.R")
source("/Users/leomarcellopoli/Documents/Bayesian/Project/BSProject2526/funzioni R/compute_metrics.R")

# ==============================================================================
# PARAMETRI E CONFIGURAZIONE
# ==============================================================================

q <- 5               # Numero di nodi del grafo
n <- 100000            # Dimensione del campione
w <- 0.6             # Prob di un arco nel DAG generato casualmente
a <- q               # Parametro (da specificare meglio l'uso)
U <- diag(1, q)      # Prior per il metodo Bayesiano
alpha_freq <- 0.01     # Soglia per il metodo PC frequentista
alpha_bayes <- 10      # Soglia per il metodo Bayesiano

# ==============================================================================
# CREAZIONE DEL DAG VERO E GENERAZIONE DATI
# ==============================================================================

# Matrice di adiacenza del DAG vero
# D0 <- matrix(c(rep(0, q),
#                rep(0, q),
#                rep(c(1, 1, 0, 0, 0), 3)), 
#              byrow = TRUE, nrow = q, ncol = q)
D0 <- rDAG(q, w)  # da BCDAG

print("DAG vero:")
D0

# Parametri della decomposizione di Cholesky modificata di Sigma
# Matrice L (lower triangular con struttura del DAG)
L <- D0 * # Matrice di adicenaza
  matrix(runif(q * q, 0.1, 1), q, q) *  # Coefficienti uniformi [0.1, 1]
  sample(c(-1, 1), size = q * q, replace = TRUE)  # Segni casuali
diag(L) <- 1

print("Matrice L:")
print(L)

# Matrice D (diagonale)
D <- diag(rep(1, q))
D

# Matrice di covarianza: Sigma = (L')^(-1) * D * L^(-1)
Sigma <- solve(t(L)) %*% D %*% solve(L)

print("Matrice di covarianza Sigma:")
print(Sigma)

# Generazione dataset multivariato normale
X <- mvtnorm::rmvnorm(n, mean = rep(0, q), sigma = Sigma)
head(X)

# ==============================================================================
# STIMA DEGLI SKELETON
# ==============================================================================

# Skeleton reale
DAG.true <- as(D0, "graphNEL")

# Metodo PC frequentista
skltn.freq <- skeleton(suffStat = list(C = cor(X), n = nrow(X)), 
                       indepTest = gaussCItest,
                       p = ncol(X),
                       alpha = alpha_freq)

# Metodo Bayesiano personalizzato
res.bayes <- my_skeleton(X = X, 
                         a = a, 
                         U = U, 
                         alpha = alpha_bayes, 
                         BayesTest = BF_Gaussian,
                         saveBF = TRUE,       
                         verbose = FALSE)
skltn.bayes <- res.bayes$object

# ==============================================================================
# VISUALIZZAZIONE RISULTATI
# ==============================================================================

#par(mfrow = c(1, 3))
#plot(DAG.true, main = "DAG Vero")
#plot(skltn.freq@graph, main = "Skeleton Frequentista")
#plot(skltn.bayes@graph, main = "Skeleton Bayesiano")
#dev.off()

D0_skeleton <- D0 + t(D0)  # Rende la matrice simmetrica
D0_skeleton[D0_skeleton > 0] <- 1  # Converte in 0/1
skltn.true <- as(D0_skeleton, "graphNEL")
#print("Shd between true skeleton and frequentist skeleton = ")
#print(shd(skltn.true, skltn.freq@graph))
#print("Shd between true skeleton and bayesian skeleton = ")
#print(shd(skltn.true, skltn.bayes@graph))

eval_freq <- evaluate_graph(skltn.true, skltn.freq@graph, undirected = TRUE)
eval_bayes <- evaluate_graph(skltn.true, skltn.bayes@graph, undirected = TRUE)
comparison <- rbind(Frequentista = eval_freq, Bayesiano = eval_bayes)
print(comparison)

#par(mfrow = c(1, 2))
#
#metrics_pr <- c("Precision", "Recall")
#colors <- c("steelblue", "coral")
#
#barplot(as.matrix(comparison[, metrics_pr]), 
#        beside = TRUE, 
#        col = colors,
#        main = "Precision e Recall",
#        ylab = "Valore",
#        ylim = c(0, 1),
#        names.arg = c("Precision", "Recall"),  
#        legend.text = rownames(comparison),
#        args.legend = list(x = "topright", cex = 0.8, bty = "n"))
#
#barplot(as.matrix(comparison[, "SHD", drop = FALSE]), 
#        beside = TRUE, 
#        col = colors,
#        main = "Structural Hamming Distance",
#        ylab = "SHD",
#        ylim = c(0, q),  # Adatta al valore massimo
#        names.arg = "SHD",  
#        legend.text = rownames(comparison),
#        args.legend = list(x = "topright", cex = 0.8, bty = "n"))  
