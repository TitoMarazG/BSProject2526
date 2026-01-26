skeleton <- function (suffStat, indepTest, alpha, labels, p, method = c("stable", 
                                                                        "original", "stable.fast"), m.max = Inf, fixedGaps = NULL, 
                      fixedEdges = NULL, NAdelete = TRUE, numCores = 1, verbose = FALSE) 
{
  cl <- match.call() # Records the function call for the output object
  
  ## --- PART 1: Input Validation & Initialization ---
  if (!missing(p)) 
    stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  
  # Ensure we have node labels; if not provided, use numbers 1 to p
  if (missing(labels)) {
    if (missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else {
    stopifnot(is.character(labels))
    if (missing(p)) p <- length(labels)
    else if (p != length(labels)) 
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
  }
  
  seq_p <- seq_len(p)
  method <- match.arg(method) # Choose between stable, original, or stable.fast
  
  # Initialize the Adjacency Matrix G (Start with a fully connected graph)
  if (is.null(fixedGaps)) {
    G <- matrix(TRUE, nrow = p, ncol = p)
  } else {
    # If fixedGaps is provided, those edges are pre-removed
    G <- !fixedGaps 
  }
  diag(G) <- FALSE # No self-loops
  
  # fixedEdges: Edges that cannot be removed by the algorithm
  if (any(is.null(fixedEdges))) {
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  
  ## --- PART 2: High-Performance Execution (C++ Call) ---
  # If stable.fast is chosen, the logic is handled by compiled C++ code
  if (method == "stable.fast") {
    indepTestName <- if (identical(indepTest, gaussCItest)) "gauss" else "rfun"
    options <- list(verbose = as.integer(verbose), 
                    m.max = as.integer(ifelse(is.infinite(m.max), p, m.max)), 
                    NAdelete = NAdelete, numCores = numCores)
    
    res <- .Call(estimateSkeleton, G, suffStat, indepTestName, 
                 indepTest, alpha, fixedEdges, options)
    
    # Reformat C++ results back into R structures
    G <- res$amat
    sepset <- lapply(seq_p, function(i) c(lapply(res$sepset[[i]], 
                                                 function(v) if (identical(v, as.integer(-1))) NULL else v), 
                                          vector("list", p - length(res$sepset[[i]]))))
    pMax <- res$pMax
    n.edgetests <- res$n.edgetests
    ord <- length(n.edgetests) - 1L
    
  } else {
    ## --- PART 3: Standard R Implementation (Iterative Testing) ---
    pval <- NULL
    sepset <- lapply(seq_p, function(.) vector("list", p)) # To store separation sets
    pMax <- matrix(-Inf, nrow = p, ncol = p) # To store maximum p-values found
    diag(pMax) <- 1
    done <- FALSE
    ord <- 0L # 'ord' is the size of the conditioning set (starting at 0)
    n.edgetests <- numeric(1)
    
    # Loop while there are still edges to test and we haven't reached max set size
    while (!done && any(G) && ord <= m.max) {
      n.edgetests[ord1 <- ord + 1L] <- 0
      done <- TRUE
      
      # Find current edges in the graph
      ind <- which(G, arr.ind = TRUE)
      ind <- ind[order(ind[, 1]), ]
      remEdges <- nrow(ind)
      
      # 'stable' method: compute neighbors once at the start of each 'ord' 
      # to avoid the removal of one edge affecting the removal of others in the same level
      if (method == "stable") {
        G.l <- split(G, gl(p, p))
      }
      
      for (i in 1:remEdges) {
        x <- ind[i, 1]
        y <- ind[i, 2]
        
        # Only test if edge exists and isn't fixed
        if (G[y, x] && !fixedEdges[y, x]) {
          nbrsBool <- if (method == "stable") G.l[[x]] else G[, x]
          nbrsBool[y] <- FALSE
          nbrs <- seq_p[nbrsBool] # Possible conditioning variables
          length_nbrs <- length(nbrs)
          
          if (length_nbrs >= ord) {
            if (length_nbrs > ord) done <- FALSE # More tests possible in next ord
            
            S <- seq_len(ord) # Initial subset of neighbors
            repeat {
              n.edgetests[ord1] <- n.edgetests[ord1] + 1
              
              # PERFORM INDEPENDENCE TEST: Are x and y independent given subset S?
              pval <- indepTest(x, y, nbrs[S], suffStat)
              
              if (is.na(pval)) pval <- as.numeric(NAdelete)
              if (pMax[x, y] < pval) pMax[x, y] <- pval
              
              # If independent (p-value > alpha), remove the edge
              if (pval >= alpha) {
                G[x, y] <- G[y, x] <- FALSE
                sepset[[x]][[y]] <- nbrs[S] # Record the set that separated them
                break
              } else {
                # Try the next combination of neighbors
                nextSet <- getNextSet(length_nbrs, ord, S)
                if (nextSet$wasLast) break
                S <- nextSet$nextSet
              }
            }
          }
        }
      }
      ord <- ord + 1L # Increment conditioning set size (0 -> 1 -> 2...)
    }
    
    # Symmetrize the pMax matrix
    for (i in 1:(p - 1)) {
      for (j in 2:p) pMax[i, j] <- pMax[j, i] <- max(pMax[i, j], pMax[j, i])
    }
  }
  
  ## --- PART 4: Output Generation ---
  # Convert result to a graph object (graphNEL)
  Gobject <- if (sum(G) == 0) {
    new("graphNEL", nodes = labels)
  } else {
    colnames(G) <- rownames(G) <- labels
    as(G, "graphNEL")
  }
  
  # Wrap everything in a pcAlgo object
  new("pcAlgo", graph = Gobject, call = cl, n = integer(0), 
      max.ord = as.integer(ord - 1), n.edgetests = n.edgetests, 
      sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
}