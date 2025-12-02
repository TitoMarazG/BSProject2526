#include "BayesCore.h"



// Implementation of Multivariate Gamma in log-scale
double log_mgamma(double x, int p) {
    // Term 1: p*(p-1)/4 * log(pi)
    double term1 = (p * (p - 1) / 4.0) * LOG_PI;

    // Term 2: Sum of lgamma terms
    double term2 = 0.0;
    for (int i = 1; i <= p; ++i) {
        // gamma(x + (1 - i)/2)
        term2 += std::lgamma(x + (1.0 - i) / 2.0);
    }

    return term1 + term2;
}



// Implementation of the normalization constant c(a, U)
// Reference: 'norm_cost' function in Function_p(XA).Rmd
double log_norm_cost(double df, const arma::mat& U) {
    int q = U.n_cols; // Dimension of the submatrix

    // Calculate log determinant of U securely using Cholesky
    // num = |U|^(a/2)  =>  log_num = (a/2) * log(|U|)
    arma::mat L;
    double log_det_U = 0.0;

    if(arma::chol(L, U)) {
        // log(det(U)) = 2 * sum(log(diag(L)))
        log_det_U = 2.0 * arma::sum(arma::log(L.diag()));
    } else {
        // Fallback for non-pd matrices (should not happen with proper priors)
        return -1e100;
    }

    double log_num = (df / 2.0) * log_det_U;

    // den = 2^(a*q/2) * Gamma_q(a/2)
    // log_den = (a*q/2)*log(2) + log_mgamma(a/2, q)
    double log_den = (df * q / 2.0) * std::log(2.0) + log_mgamma(df / 2.0, q);

    // Result = num / den => log_num - log_den
    return log_num - log_den;
}



// Implementation of p(XA)
// Reference: 'p_XA' function in Function_p(XA).Rmd
double log_p_XA(const arma::mat& S, int n, double a, const arma::mat& U, const arma::uvec& A) {
    // Handle empty set case (A is empty) -> standard convention p(X_empty) = 1 -> log = 0
    if (A.n_elem == 0) return 0.0;

    int q_total = S.n_cols;
    int card_A = A.n_elem;
    int card_nonA = q_total - card_A; // Corresponds to 'card_nonA' in R

    // Extract submatrices
    // S_AA = S[A, A]
    // U_AA = U[A, A]
    arma::mat S_AA = S.submat(A, A);
    arma::mat U_AA = U.submat(A, A);
    arma::mat Post_AA = U_AA + S_AA; // U + S

    // Calculate priors and posteriors normalization constants
    // R code: c_prior <- norm_cost(a - card_nonA, U_AA)
    double log_c_prior = log_norm_cost(a - card_nonA, U_AA);

    // R code: c_post <- norm_cost(a + n - card_nonA, U_AA + S_AA)
    double log_c_post  = log_norm_cost(a + n - card_nonA, Post_AA);

    // Calculate factor
    // R code: factor = (2 * pi)^(- (n * card_A) / 2)
    double log_factor = - (n * card_A / 2.0) * std::log(2 * arma::datum::pi);

    // Final result: factor * (c_prior / c_post)
    return log_factor + log_c_prior - log_c_post;
}



double log_p_DAG(const arma::mat& S, int n, const arma::mat& adj_matrix, double a, const arma::mat& U) {
    int n_nodes = adj_matrix.n_cols;
    double total_score = 0.0;

    for (int i = 0; i < n_nodes; ++i) {
        // Find parents: indices where column i has value 1
        arma::uvec pa = arma::find(adj_matrix.col(i) == 1);

        // Family: parents + current node i
        arma::uvec fa(pa.n_elem + 1);
        fa.head(pa.n_elem) = pa;
        fa(pa.n_elem) = i;

        // Score component: log P(X_fa) - log P(X_pa)
        double score_fa = log_p_XA(S, n, a, U, fa);
        double score_pa = 0.0;

        if (pa.n_elem > 0) {
            score_pa = log_p_XA(S, n, a, U, pa);
        }

        total_score += (score_fa - score_pa);
    }
    return total_score;
}



double calculate_BF(const arma::mat& S, int n, const arma::mat& adj0, const arma::mat& adj1,
                    double a, const arma::mat& U, int node_idx) {
    // Extract families for the specific node in both graphs

    // Graph 0
    arma::uvec pa0 = arma::find(adj0.col(node_idx) == 1);
    arma::uvec fa0(pa0.n_elem + 1); fa0.head(pa0.n_elem) = pa0; fa0(pa0.n_elem) = node_idx;
    double log_local_score0 = log_p_XA(S, n, a, U, fa0) - (pa0.n_elem > 0 ? log_p_XA(S, n, a, U, pa0) : 0.0);

    // Graph 1
    arma::uvec pa1 = arma::find(adj1.col(node_idx) == 1);
    arma::uvec fa1(pa1.n_elem + 1); fa1.head(pa1.n_elem) = pa1; fa1(pa1.n_elem) = node_idx;
    double log_local_score1 = log_p_XA(S, n, a, U, fa1) - (pa1.n_elem > 0 ? log_p_XA(S, n, a, U, pa1) : 0.0);

    // Bayes Factor = P(D0) / P(D1)
    // log(BF) = log(P0) - log(P1)
    return std::exp(log_local_score0 - log_local_score1);
}