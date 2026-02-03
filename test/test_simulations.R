# SETUP INIZIALE================================================================

setwd("~/Documents/Bayesian/Project/BSProject2526")

library(pcalg)
library(BiocGenerics)
library(BCDAG)

rm(list = ls()) 
graphics.off() 

source("/Users/leomarcellopoli/Documents/Bayesian/Project/BSProject2526/funzioni R/personalized_my_skeleton_senza_stable_fast.R")
source("/Users/leomarcellopoli/Documents/Bayesian/Project/BSProject2526/funzioni R/compute_metrics.R")

# Liste di parametri da testare
q_list <- c(5, 10, 20)
n_list <- c(100, 500, 1000)
repetitions <- 1:20
w <- 0.6
alpha_freq <- 0.01
alpha_bayes <- 10

# Crea lista vuota per salvare i risultati
results_list <- list()
counter <- 1

# LOOP ESTERNO: varia q
for (q in q_list) {
  
  cat("\n========================================\n")
  cat("Testing with q =", q, "nodes\n")
  cat("========================================\n")
  
  # LOOP INTERNO: varia n
  for (n in n_list) {
    
    # LOOP repetitions (stessa (n, q) per #repetitions volte)
    
    for (rep in repetitions) {
      
      cat("  Running simulation: q =", q, ", n =", n, ", repetition =", rep, "... ")
      
      # Genera DAG e dati
      set.seed(1 + rep)  # Seed diverso per ogni simulazione
      D0 <- rDAG(q, w)
      U <- diag(1, q)
      
      L <- D0 * matrix(runif(q * q, 0.1, 1), q, q) * 
        sample(c(-1, 1), size = q * q, replace = TRUE)
      diag(L) <- 1
      D <- diag(rep(1, q))
      Sigma <- solve(t(L)) %*% D %*% solve(L)
      
      X <- mvtnorm::rmvnorm(n, mean = rep(0, q), sigma = Sigma)
      
      # Skeleton vero
      D0_skeleton <- D0 + t(D0)
      D0_skeleton[D0_skeleton > 0] <- 1
      skltn.true <- as(D0_skeleton, "graphNEL")
      
      # Stima frequentista
      skltn.freq <- skeleton(suffStat = list(C = cor(X), n = nrow(X)), 
                             indepTest = gaussCItest,
                             p = ncol(X),
                             alpha = alpha_freq)
      
      # Stima bayesiana
      #res.bayes <- my_skeleton(X = X, a = q, U = U, 
      #                         alpha = alpha_bayes, 
      #                         BayesTest = BF_Gaussian,
      #                         saveBF = TRUE, verbose = FALSE)
      #skltn.bayes <- res.bayes$object
      skltn.bayes <- my_skeleton(X = X, a = q, U = U, 
                                 alpha = alpha_bayes, 
                                 BayesTest = BF_Gaussian,
                                 saveBF = FALSE, verbose = FALSE)
      
      # Calcola metriche
      eval_freq <- evaluate_graph(skltn.true, skltn.freq@graph, undirected = TRUE)
      eval_bayes <- evaluate_graph(skltn.true, skltn.bayes@graph, undirected = TRUE)
      
      # Salva risultati
      results_list[[counter]] <- data.frame(
        q = q,
        n = n,
        method = c("Frequentist", "Bayesian"),
        SHD = c(eval_freq$SHD, eval_bayes$SHD),
        Precision = c(eval_freq$Precision, eval_bayes$Precision),
        Recall = c(eval_freq$Recall, eval_bayes$Recall),
        F1_Score = c(eval_freq$F1_Score, eval_bayes$F1_Score),
        TP = c(eval_freq$TP, eval_bayes$TP),
        FP = c(eval_freq$FP, eval_bayes$FP),
        FN = c(eval_freq$FN, eval_bayes$FN),
        TN = c(eval_freq$TN, eval_bayes$TN),
        test_dataset = repetitions[rep]
      )
      
      cat("Done!\n")
      counter <- counter + repetitions[rep]
    }
  }
}

# Risultati in un unico data frame
results_df <- do.call(rbind, results_list)
rownames(results_df) <- NULL

View(results_df)

#write.csv(results_df, "simulation_results.csv", row.names = FALSE)