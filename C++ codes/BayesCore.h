#ifndef BAYESCORE_H
#define BAYESCORE_H

#include <armadillo>
#include <cmath>



// Constant for Log(Pi)
const double LOG_PI = std::log(arma::datum::pi);



/**
 * Calculates the Log Multivariate Gamma function.
 * Represents the log of the term 'mvgamma_direct' in R.
 * - @param x The argument (usually degrees of freedom / 2)
 * - @param p The dimension
 * - @return Log-Gamma value
 */
double log_mgamma(double x, int p);



/**
 * Calculates the log of the normalization constant c(a, U).
 * Implements the formula from 'norm_cost' in Function_p(XA).Rmd using log-scale.
 * Formula: log(|U|^(a/2)) - log(2^(a*p/2) * Gamma_p(a/2))
 * - @param df Degrees of freedom (parameter 'a')
 * - @param U Scale matrix
 * - @return Log normalization constant
 */
double log_norm_cost(double df, const arma::mat& U);



/**
 * Calculates the Marginal Likelihood p(XA) in log-scale.
 * Corresponds to 'p_XA' in your R script.
 * - @param S Scatter matrix X'X (pre-computed)
 * - @param n Number of observations
 * - @param a Hyperparameter (degrees of freedom)
 * - @param U Hyperparameter (scale matrix)
 * - @param A Vector of variable indices (0-based)
 * - @return Log p(XA)
 */
double log_p_XA(const arma::mat& S, int n, double a, const arma::mat& U, const arma::uvec& A);
/* Reference on R: p_XA <- function(X, a, U, A) */



/**
 * Calculates the total log-score of a DAG (BGe Score).
 * Iterates over nodes to compute sum(log_p_XA(fa) - log_p_XA(pa)).
 */
double log_p_DAG(const arma::mat& S, int n, const arma::mat& adj_matrix, double a, const arma::mat& U);
/* Reference on R: P_DAG <- function(X, DAG, a, U) */



/**
 * Calculates the Bayes Factor between two DAGs for a specific node.
 * Returns exp(Score0 - Score1).
 */
double calculate_BF(const arma::mat& S, int n, const arma::mat& adj0, const arma::mat& adj1,
                    double a, const arma::mat& U, int node_idx);
/* Reference on R: BF <- function(S, DAGM0, DAGM1, a, U, nodo_da_testare) */

#endif // BAYESCORE_H