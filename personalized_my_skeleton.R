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
my_skeleton <- function (X, a, U,  # QUESTI SONO GLI INPUT PER IL CASO BAYESIANO
                         BayesTest = BF_Gaussian,  # CHE PRIMA ERA indepTest
                         alpha,  # DIVENTA LA SOGLIA
                         #labels, p,  # QUESTI NON SERVONO PIù DATO CHE DIAMO IN INPUT X (cioè XX e n)
                         method = c("stable", "original", "stable.fast"),  # stable.fast servirà per quando useremo C++
                         # per ora lavoriamo su stable, cioè solo in R
                         m.max = Inf,
                         fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, numCores = 1, verbose = FALSE) 
{
  cl <- match.call() # serve a catturare e memorizzare la "chiamata" esatta che è stata fatta alla funzione corrente,
  # prende il modo in cui l'utente ha scritto la funzione (inclusi gli argomenti) e lo trasforma 
  # in un oggetto manipolabile all'interno del codice.
  
  # CHECK INPUT
  if (missing(X))
    stop("La matrice X è obbligatoria!")
  if (is.data.frame(X))
    X <- as.matrix(X)
  
  # CALCOLO STATISTICHE SUFFICIENTI
  XX <- t(X) %*% X
  n <- nrow(X)
  p <- ncol(X)
  labels <- colnames(X)
  if (is.null(labels))
    labels <- as.character(seq_len(p))
  
  seq_p <- seq_len(p) # vettore che va da 1 a p
  method <- match.arg(method) # serve a controllare che l'utente abbia inserito un'opzione valida
  # e a completarla automaticamente se l'ha scritta solo in parte,
  # oppure a scegliere la predefinita (la prima).
  
  # DA QUI ...
  if (is.null(fixedGaps))
    G <- matrix(TRUE, nrow = p, ncol = p)
  else if (!identical(dim(fixedGaps), c(p, p))) 
    stop("Dimensions of the dataset and fixedGaps do not agree.")
  else if (!identical(fixedGaps, t(fixedGaps))) 
    stop("fixedGaps must be symmetric")
  else
    G <- !fixedGaps
  diag(G) <- FALSE
  
  if (any(is.null(fixedEdges))) {
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedEdges), c(p, p))) 
    stop("Dimensions of the dataset and fixedEdges do not agree.")
  else if (!identical(fixedEdges, t(fixedEdges))) 
    stop("fixedEdges must be symmetric")
  # ... A QUI: inizializza la matrice G (che rappresenta il grafo completo)
  
  stopifnot((is.integer(numCores) || is.numeric(numCores)) && numCores > 0)
  if (numCores > 1 && method != "stable.fast") {
    warning("Argument numCores ignored: parallelization only available for method = 'stable.fast'")
  }
  
  # ## STABLE FAST ## CIOÈ LA VERSIONE CHE RICHIAMA LA FUNZIONE DA C++
  if (method == "stable.fast") {
    if (identical(indepTest, gaussCItest)) 
      indepTestName <- "gauss"
    else indepTestName <- "rfun"
    options <- list(verbose = as.integer(verbose), m.max = as.integer(ifelse(is.infinite(m.max), 
                                                                             p, m.max)), NAdelete = NAdelete, numCores = numCores)
    
    # QUI VIENE CHIAMATA LA FUNZIONE CHE IMPLEMENTA L'ALGORITMO PC STABLE FAST DA C++
    res <- .Call(estimateSkeleton, G, suffStat, indepTestName, 
                 indepTest, alpha, fixedEdges, options)
    
    G <- res$amat # ADJANCENCY MATRIX DEL GRAFO STIMATO DA C++
    
    # SETTO LE VARIABILI DA RITORNARE
    sepset <- lapply(seq_p, function(i) c(lapply(res$sepset[[i]], # SEPSET SONO I SEPARATION SET RICAVATI DALL'OUTPUT DELLA FUNZIONE C++ 
                                                 function(v) if (identical(v, as.integer(-1))) NULL else v), 
                                          vector("list", p - length(res$sepset[[i]]))))
    pMax <- res$pMax  # P-VALUES MASSIMI PER OGNI COPPIA DI VARIABILI RICAVATI DALL'OUTPUT DELLA FUNZIONE C++
    n.edgetests <- res$n.edgetests  # NUMERO DI TEST ESEGUITI AD OGNI ORDINE
    ord <- length(n.edgetests) - 1L  # MASSIMO ORDINE RAGGIUNTO
  }
  # QUI FINISCE stable.fast, DOVREMO IMPLEMENTARLO DA C++
  
  else {
    # ## NOT STABLE FAST ##
    # ORA RIMANGO DENTRO R: QUI INIZIA L'IMPLEMENTAZIONE DELL'ALGORITMO PC STABLE ORIGINALE E STABLE
    
    # SETTO LE VARIABILI PER FARE: IL CICLO WHILE (CUORE DELL'ALGORITMO)
    sepset <- lapply(seq_p, function(.) vector("list", p))  # SEPARATION SETS
    valMax <- matrix(-Inf, nrow = p, ncol = p)  # MATRICE DEI P-VALUES MASSIMI
    diag(valMax) <- Inf  # LE DIAGONALI SONO 1 (MA VANNO RAGIONATI)
    
    done <- FALSE  # VARIABILE BOOLEANA PER USCIRE DAL CICLO WHILE
    ord <- 0L  # ORD SAREBBE LA DIMENSIONE DEI CONDITIONING SETS
    n.edgetests <- numeric(1)  # VETTORE CHE CONTIENE IL NUMERO DI TEST ESEGUITI AD OGNI ORDINE
    # FINE SETTING
    
    # CUORE DELL'ALGORITMO: CICLO WHILE 
    while (!done && any(G) && ord <= m.max) {  # ORD SAREBBE LA DIMENSIONE DEI CONDITIONING SETS
      
      n.edgetests[ord1 <- ord + 1L] <- 0 
      done <- TRUE
      
      # QUI IDENTIFICO GLI ARCHI ATTUALI
      ind <- which(G, arr.ind = TRUE)
      ind <- ind[order(ind[, 1]), ]
      remEdges <- nrow(ind)
      
      #VERBOSE
      if (verbose) 
        cat("Order=", ord, "; remaining edges:", remEdges, 
            "\n", sep = "")
      
      #METHOD STABLE, VA INDAGATO O CHIEDERE AL PROF SE RIMANERE SU STABLE O ORIGINAL
      if (method == "stable") {
        G.l <- split(G, gl(p, p)) 
      }
      
      for (i in 1:remEdges) {  # SCORRIAMO I REMAINING EDGES
        
        #VERBOSE
        if (verbose && (verbose >= 2 || i%%100 == 0)) 
          cat("|i=", i, "|iMax=", remEdges, "\n")
        
        x <- ind[i, 1]
        y <- ind[i, 2]
        if (G[y, x] && !fixedEdges[y, x]) {
          
          # CREAZIONE DEI VICINI nbrsBool E DI length_nbrs
          nbrsBool <- if (method == "stable") G.l[[x]] else G[, x]
          nbrsBool[y] <- FALSE
          nbrs <- seq_p[nbrsBool]      #vettore degli indici dei vicini
          length_nbrs <- length(nbrs)
          
          if (length_nbrs >= ord) {  #  SE FOSSE length_nbrs < ord NON SI ENTRA NEL PCALG
            if (length_nbrs > ord) done <- FALSE 
            
            S <- seq_len(ord)  # CONDITIONING SET (non sono le variabili, ma gli indici dentro nbrs)
            # sono i sottoinsiemi dei vicini di dimensione "ord"
            
            repeat {
              n.edgetests[ord1] <- n.edgetests[ord1] + 1
              
              # CHIAMATA DEL BF TEST
              BFval <- BayesTest(XX, n, a, U, x, y, nbrs[S])
              
              #VERBOSE
              if (verbose) 
                cat("x=", x, " y=", y, " S=", nbrs[S], 
                    ": BF =", BFval, "\n")
              
              if (is.na(BFval))  # QUESTO SE IL TEST FALLISCE
                BFval <- as.numeric(NAdelete)
              
              if (valMax[x, y] < BFval) 
                valMax[x, y] <- BFval  #Salva il massimo BFvalue osservato finora per quella coppia (x,y) tra i vari conditioning set provati.
              
              # CONTROLLARE COL PROF, DUBBIO DELLA LOGICA, FORSE VANNO SALVATI I valmin?
              if (BFval >= alpha) {
                G[x, y] <- G[y, x] <- FALSE   # l’arco non deve esserci nello skeleton, lo elimino in G
                sepset[[x]][[y]] <- nbrs[S]   # salvo il separating set
                break
              }
              else {
                nextSet <- getNextSet(length_nbrs, ord, S)  # genera la prossima combinazione di indici
                if (nextSet$wasLast) 
                  break
                S <- nextSet$nextSet   # aggiorni e e ripeti il test
              }
            }
          }
        }
      }
      ord <- ord + 1L  # RIFACCIO TUTTO CON L'ORDINE DEI SEPARATION AUMENTATO DI 1
    }
    # QUI RENDO SIMMETRICA LA MATRICE DEI BF MASSIMI, DUBBIO?
    for (i in 1:(p - 1)) {
      for (j in 2:p) valMax[i, j] <- valMax[j, i] <- max(valMax[i, j], valMax[j, i])
    }
  }
  
  # QUI SI COSTRUSCE L'OGGETTO GRAFO E RITORNA IL RISULTATO
  Gobject <- if (sum(G) == 0) {
    new("graphNEL", nodes = labels)
  }
  else {
    colnames(G) <- rownames(G) <- labels
    as(G, "graphNEL")
  }
  new("pcAlgo", graph = Gobject, call = cl, n = integer(0), 
      max.ord = as.integer(ord - 1), n.edgetests = n.edgetests, 
      sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
}

