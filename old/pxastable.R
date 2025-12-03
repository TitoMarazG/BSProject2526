library(mvtnorm)

# multivariate log-gamma
lmgamma_p <- function(a, p) (p*(p-1)/4)*log(pi) + sum(lgamma(a + (1 - 1:p)/2))

# stable log-determinant via Cholesky
logdet_chol <- function(M) 2 * sum(log(diag(chol(M))))

# synthetic, numerically stable p(X_A)
p_XA <- function(X, a, U, A) {
  n <- nrow(X)
  q <- ncol(X)
  A <- sort(as.integer(A))
  p <- length(A)
  card_nonA <- q - p
  
  y <- (a - card_nonA)/2
  z <- y + n/2
  if (y <= (p-1)/2) stop("Invalid hyperparameters: y too small for multivariate gamma")
  
  S <- t(X) %*% X
  U_AA <- U[A,A,drop=FALSE]
  S_AA <- S[A,A,drop=FALSE]
  
  # coupled log-space terms
  paired1 <- y * logdet_chol(U_AA) - lmgamma_p(y, p)
  paired2 <- -z * logdet_chol(U_AA + S_AA) + lmgamma_p(y + n/2, p)
  
  logp <- -(n*p)/2 * log(pi) + paired1 + paired2
  list(logp = as.numeric(logp), p = exp(logp))
}

# --- Example usage ---
set.seed(99)
q <- 5; n <- 200; a <- q; A <- 1:3
U <- diag(1, q)
L <- matrix(runif(q*q,1,3),q,q)*sample(c(-1,1),q*q,TRUE)
diag(L) <- 1
D <- diag(1,q)
Sigma <- solve(t(L)) %*% D %*% solve(L)
X <- rmvnorm(n, rep(0,q), Sigma)

res <- p_XA(X, a, U, A)
res$logp  # log p(X_A)
res$p     # p(X_A)
