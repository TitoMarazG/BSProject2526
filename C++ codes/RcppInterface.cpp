#include <RcppArmadillo.h>
#include "BayesCore.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Wrapper for p_XA (returns log value)
// [[Rcpp::export]]
double rcpp_log_p_XA(arma::mat S, int n, double a, arma::mat U, arma::uvec A) {
    // Convert R indices (1-based) to C++ indices (0-based)
    A = A - 1;
    return log_p_XA(S, n, a, U, A);
}

// Wrapper for full DAG score
// [[Rcpp::export]]
double rcpp_log_p_DAG(arma::mat S, int n, arma::mat adj_matrix, double a, arma::mat U) {
    return log_p_DAG(S, n, adj_matrix, a, U);
}

// Wrapper for Bayes Factor
// [[Rcpp::export]]
double rcpp_BF(arma::mat S, int n, arma::mat adj0, arma::mat adj1, double a, arma::mat U, int node_idx) {
    // Convert node index to 0-based
    return calculate_BF(S, n, adj0, adj1, a, U, node_idx - 1);
}