# ==============================================================================
# 1. FUNZIONI AUSILIARIE (BF e calcolo BGe Score)
# ==============================================================================
# Logica: BF > soglia implica che u non aggiunge info a v dato S -> Indipendenza
BF_Gaussian <- function(XX, n, a, U, u, v, S){
  return (exp(log.P_DAGM_puntuale(XX, n, a, U, v, S) - log.P_DAGM_puntuale(XX,n, a, U, v, c(S, u))))
}

log.P_DAGM_puntuale <- function(XX, n, a, U, v, S) {
  return(log.p_XA(XX, n, a, U, c(S, v)) - log.p_XA(XX, n, a, U, S))
}

log.p_XA <- function(XX, n, a, U, A){
  nonA <- setdiff(1:ncol(XX), A)
  A.c <- length(A)  # cardinalità di A
  nonA.c <- length(nonA)  # cardinalità di nonA
  
  # Gestione indici: U e XX devono essere subsettabili
  UAA <- U[A, A, drop = FALSE] 
  XXAA <- XX[A, A, drop = FALSE]
  
  return(log(2*pi) * (-n*A.c/2)
         + log.norm_cost(a-nonA.c, UAA)
         - log.norm_cost(a+n-nonA.c, UAA+XXAA))
}

log.norm_cost <- function(a, U){
  q = ncol(U)
  return(a/2*log(det(U))
         - a*q/2*log(2)
         - log.multivariate_gamma(q, a/2))
}

log.multivariate_gamma <- function(p, x){
  val = 0
  for(j in 1:p)
    val = val + lgamma(x+(1-j)/2)
  return(p*(p-1)/4*log(pi) + val)
}

# ==============================================================================
# 2. FUNZIONE getNextSet (Necessaria per le combinazioni)
# ==============================================================================
getNextSet <- function (n, k, set) {
  chInd <- k - (zeros <- sum((seq(n - k + 1, n) - set) == 0))
  wasLast <- (chInd == 0)
  if (!wasLast) {
    set[chInd] <- s.ch <- set[chInd] + 1
    if (chInd < k) 
      set[(chInd + 1):k] <- seq(s.ch + 1L, s.ch + zeros)
  }
  list(nextSet = set, wasLast = wasLast)
}

# ==============================================================================
# 3. MY_SKELETON (Bayesian Version)
# ==============================================================================
my_skeleton <- function (X, a, U,  # INPUT BAYESIANI
                         BayesTest = BF_Gaussian,  # Default alla tua funzione
                         alpha,  # SOGLIA BF
                         method = c("BayesTest", "original"), 
                         m.max = Inf,
                         fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, numCores = 1, verbose = FALSE) 
{
  cl <- match.call() 
  
  # CHECK INPUT
  if (missing(X)) stop("La matrice X è obbligatoria!")
  if (is.data.frame(X)) X <- as.matrix(X)
  
  # CALCOLO STATISTICHE SUFFICIENTI
  XX <- t(X) %*% X
  n <- nrow(X)
  p <- ncol(X)
  labels <- colnames(X)
  if (is.null(labels)) labels <- as.character(seq_len(p))
  
  seq_p <- seq_len(p) 
  method <- match.arg(method) 
  
  # INIZIALIZZAZIONE G
  if (is.null(fixedGaps))
    G <- matrix(TRUE, nrow = p, ncol = p)
  else {
    if (!identical(dim(fixedGaps), c(p, p))) stop("Dimensions of fixedGaps do not agree.")
    if (!identical(fixedGaps, t(fixedGaps))) stop("fixedGaps must be symmetric")
    G <- !fixedGaps
  }
  diag(G) <- FALSE
  
  if (any(is.null(fixedEdges))) {
    fixedEdges <- matrix(FALSE, nrow = p, ncol = p)
  } else {
    if (!identical(dim(fixedEdges), c(p, p))) stop("Dimensions of fixedEdges do not agree.")
    if (!identical(fixedEdges, t(fixedEdges))) stop("fixedEdges must be symmetric")
  }
  
  stopifnot((is.integer(numCores) || is.numeric(numCores)) && numCores > 0)
  
  # -------------------------------------------------------------------------
  # BLOCCO STABLE.FAST (Disattivato correttamente con commenti)
  # -------------------------------------------------------------------------
  # if (method == "stable.fast") {
  #   stop("Il metodo stable.fast richiede implementazione C++ non disponibile.")
  # }
  # -------------------------------------------------------------------------
  
  # BLOCCO R PURO (Stable / Original)
  
  # Setup variabili ciclo
  sepset <- lapply(seq_p, function(.) vector("list", p)) 
  
  # valMax: tiene traccia della MASSIMA evidenza di indipendenza trovata.
  valMax <- matrix(-Inf, nrow = p, ncol = p) 
  diag(valMax) <- 1 # Diagonale ininfluente
  
  done <- FALSE 
  ord <- 0L 
  n.edgetests <- numeric(1) 
  
  # CUORE DELL'ALGORITMO: CICLO WHILE 
  while (!done && any(G) && ord <= m.max) { 
    
    n.edgetests[ord1 <- ord + 1L] <- 0 
    done <- TRUE
    
    ind <- which(G, arr.ind = TRUE)
    ind <- ind[order(ind[, 1]), ]
    remEdges <- nrow(ind)
    
    if (verbose) cat("Order=", ord, "; remaining edges:", remEdges, "\n", sep = "")
    
    # METHOD STABLE: Snapshot dei vicini
    if (method == "stable") {
      G.l <- split(G, gl(p, p)) 
    }
    
    for (i in 1:remEdges) { 
      
      if (verbose && (verbose >= 2 || i%%100 == 0)) cat("|i=", i, "|iMax=", remEdges, "\n")
      
      x <- ind[i, 1]
      y <- ind[i, 2]
      
      if (G[y, x] && !fixedEdges[y, x]) {
        
        nbrsBool <- if (method == "stable") G.l[[x]] else G[, x]
        nbrsBool[y] <- FALSE
        nbrs <- seq_p[nbrsBool]      
        length_nbrs <- length(nbrs)
        
        if (length_nbrs >= ord) { 
          if (length_nbrs > ord) done <- FALSE 
          
          S <- seq_len(ord) 
          
          repeat {
            n.edgetests[ord1] <- n.edgetests[ord1] + 1
            
            # CHIAMATA DEL BF TEST (Corretto: usa BayesTest e i nuovi input)
            BFval <- BayesTest(XX, n, a, U, x, y, nbrs[S])
            
            if (verbose) cat("x=", x, " y=", y, " S=", nbrs[S], ": BF =", BFval, "\n")
            
            if (is.na(BFval)) BFval <- as.numeric(NAdelete)
            
            # Salviamo il massimo BF (massima evidenza indipendenza)
            if (valMax[x, y] < BFval) valMax[x, y] <- BFval 
            
            # DECISIONE: Se BF > alpha (Soglia), assumiamo indipendenza -> Taglio arco
            if (BFval >= alpha) {
              G[x, y] <- G[y, x] <- FALSE 
              sepset[[x]][[y]] <- nbrs[S] 
              break
            } else {
              nextSet <- getNextSet(length_nbrs, ord, S) 
              if (nextSet$wasLast) break
              S <- nextSet$nextSet 
            }
          }
        }
      }
    }
    ord <- ord + 1L 
  }
  
  # Simmetrizzazione valMax
  for (i in 1:(p - 1)) {
    for (j in 2:p) valMax[i, j] <- valMax[j, i] <- max(valMax[i, j], valMax[j, i])
  }
  
  # OUTPUT PACKAGING
  Gobject <- if (sum(G) == 0) {
    new("graphNEL", nodes = labels)
  } else {
    colnames(G) <- rownames(G) <- labels
    as(G, "graphNEL") 
  }
  
  # Creazione oggetto pcAlgo.
  # pMax = valMax (Corretto: passiamo i BF massimi nello slot pMax)
  new("pcAlgo", graph = Gobject, call = cl, n = integer(0), 
      max.ord = as.integer(ord - 1), n.edgetests = n.edgetests, 
      sepset = sepset, pMax = valMax, zMin = matrix(NA, 1, 1))
}