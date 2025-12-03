#include <iostream>
#include <armadillo>
#include "BayesCore.h"

int main() {
    // 1. SIMULATION PARAMETERS (Replica of the R script)
    int q = 5;      // Number of nodes
    int n = 200;    // Number of simulations

    // Set seed for reproducibility
    arma::arma_rng::set_seed(42);

    std::cout << "--- Start Data Simulation ---" << std::endl;

    // 2. CONSTRUCTION OF TRUE DAG (D0)
    // R: matrix(c(rep(0,q), rep(0, q), rep(c(1,1,0,0,0),3)), byrow = T ...)
    // Row 1 (idx 0): 0 0 0 0 0
    // Row 2 (idx 1): 0 0 0 0 0
    // Row 3 (idx 2): 1 1 0 0 0
    // Row 4 (idx 3): 1 1 0 0 0
    // Row 5 (idx 4): 1 1 0 0 0
    arma::mat D0(q, q, arma::fill::zeros);
    arma::rowvec parents_pattern = {1, 1, 0, 0, 0};
    D0.row(2) = parents_pattern;
    D0.row(3) = parents_pattern;
    D0.row(4) = parents_pattern;

    std::cout << "True DAG (D0):" << std::endl;
    D0.print();

    // 3. PARAMETER GENERATION (L and Sigma)
    // R: L = D0 * runif(1,3) * sample(-1,1); diag(L) = 1
    arma::mat L(q, q, arma::fill::eye); // Diagonal set to 1 by default

    for(int i=0; i<q; ++i) {
        for(int j=0; j<q; ++j) {
            if(D0(i,j) == 1) {
                // Generate random weight between 1 and 3
                double weight = (arma::randu() * 2.0) + 1.0;
                // Generate sign (-1 or 1)
                double sign = (arma::randu() > 0.5) ? 1.0 : -1.0;
                L(i,j) = weight * sign;
            }
        }
    }

    // Calculate Sigma = (L')^-1 * D * (L)^-1
    // Since D is identity, Sigma = inv(L.t()) * inv(L)
    arma::mat L_inv = arma::inv(L);
    arma::mat Sigma = L_inv.t() * L_inv;

    // 4. DATA SIMULATION X
    // X ~ N(0, Sigma)
    arma::vec mu(q, arma::fill::zeros);
    // mvnrnd generates (q x n), so we transpose to get (n x q) as in R
    arma::mat X = arma::mvnrnd(mu, Sigma, n).t();

    // Calculate S (Scatter Matrix X'X) required by your functions
    arma::mat S = X.t() * X;

    std::cout << "Generated data (first 5 rows):" << std::endl;
    X.head_rows(5).print();
    std::cout << "Scatter Matrix S (X'X) calculated." << std::endl;


    // 5. C** FUNCTIONS TESTING
    std::cout << "\n--- BayesCore Function Testing ---" << std::endl;

    // Hyperparameters (Priors)
    double a = q + 2;
    arma::mat U = arma::eye(q, q);
    arma::mat D_empty(q, q, arma::fill::zeros);

    // --- GLOBAL SCORE TEST ---
    double score_D0 = log_p_DAG(S, n, D0, a, U);
    double score_empty = log_p_DAG(S, n, D_empty, a, U);

    std::cout << "Global Log-Score (D0):    " << score_D0 << std::endl;
    std::cout << "Global Log-Score (Empty): " << score_empty << std::endl;


    // --- NODE-WISE BAYES FACTOR LOOP ---
    std::cout << "\n--- Node-wise Bayes Factor Analysis (D0 vs Empty) ---" << std::endl;
    std::cout << "Comparing True DAG (Model 0) vs Empty DAG (Model 1) for each node." << std::endl;
    std::cout << "BF > 1 favors True DAG. BF < 1 favors Empty DAG." << std::endl;
    std::cout << "----------------------------------------------------------------" << std::endl;

    for (int i = 0; i < q; ++i) {
        // Calculate Bayes Factor for node i
        double bf = calculate_BF(S, n, D0, D_empty, a, U, i);

        // Identify true parents in D0 for display
        arma::uvec parents = arma::find(D0.col(i) == 1);

        std::cout << "Node " << i << " | True Parents: { ";
        if(parents.n_elem == 0) std::cout << "- ";
        else for(auto p : parents) std::cout << p << " ";
        std::cout << "} ";

        // Output BF with formatting
        std::cout << "| BF: " << std::setw(10) << bf;

        // Interpret result
        if (parents.n_elem == 0) {
            // For root nodes, D0 and Empty are structurally identical (0 parents).
            // BF should be exactly 1.0 (or extremely close due to floating point).
            if (std::abs(bf - 1.0) < 1e-9) {
                std::cout << "  [OK] (Models are identical)";
            } else {
                std::cout << "  [?] (Unexpected: should be 1.0)";
            }
        } else {
            // For children nodes, D0 has parents, Empty does not.
            // We expect BF > 1 if the data supports the dependencies.
            if (bf > 1.0) {
                std::cout << "  [OK] (Evidence for Parents)";
            } else {
                std::cout << "  [WARNING] (Evidence against Parents)";
            }
        }
        std::cout << std::endl;
    }
    std::cout << "----------------------------------------------------------------" << std::endl;
/*
    // Bayes Factor for a specific node
    // Testing node 2 (which has parents 0 and 1 in D0)
    // Compare D0 against D_empty for node 2
    int node_idx = 2; // 0-based index (corresponds to 3rd node in R)
    double bf = calculate_BF(S, n, D0, D_empty, a, U, node_idx);

    std::cout << "\nBayes Factor for node " << node_idx << " (D0 vs Empty): " << bf << std::endl;
    if (bf > 1.0) {
        std::cout << ">> OK: BF > 1, evidence in favor of the correct parents." << std::endl;
    } else {
        std::cout
*/
    return 0;
}