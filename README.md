# Bayesian Constraint-Based Structure Learning

Project for the **Bayesian Statistics** course, A.Y. 2025/2026.

## Overview

This project implements a **Bayesian alternative to the PC algorithm** (Spirtes et al., 2000) for learning the structure of Directed Acyclic Graphs (DAGs) from data. Standard constraint-based structure learning relies on frequentist hypothesis tests (e.g. Fisher's z-test), which depend on asymptotic approximations, give no uncertainty quantification, and tend to produce over-connected graphs in high-dimensional settings (q ≫ n).

We replace the frequentist conditional independence tests with **Bayes Factors**, computed in closed form under a conjugate Gaussian–Wishart model, giving:

- A fully Bayesian constraint-based structure learning algorithm
- Posterior probabilities of edge inclusion, instead of binary accept/reject decisions
- Better empirical performance (lower Structural Hamming Distance) in high-dimensional regimes

## Method

### Model

Data are assumed multivariate Gaussian, `x_1, ..., x_n | Ω ~ N_q(0, Ω⁻¹)`, with a conjugate Wishart prior `Ω ~ W(a, U)` on the precision matrix. Under a DAG structure `D`, the marginal likelihood factorizes over node families:

```
p(X | D) = Π_j  p(X_fa(j)) / p(X_pa(j))
```

where `pa(j)` are the parents of node `j` and `fa(j) = {j} ∪ pa(j)`. Each local marginal `p(X_A)` has a closed form via the Wishart normalizing constant, making the score fast to compute for any subset.

### Bayesian conditional independence test

For an edge `u — v` given a conditioning set `S`, we compare `H0: u ⊥ v | S` against `H1: u ⊥̸ v | S` via the Bayes Factor:

```
BF01 = p(X | D0) / p(X | D1)
```

computed locally from the closed-form marginals above, avoiding a full graph rescore for every candidate edge.

### Algorithm: Stable Bayesian PC

An order-independent ("stable") version of the PC algorithm:

1. Initialize a complete undirected graph.
2. For increasing conditioning-set size `l = 0, 1, 2, ...`, test each adjacent pair `(u, v)` against all subsets `S ⊆ adj(u) \ {v}` of size `l` using the Bayes Factor.
3. Remove edge `(u, v)` and record the separator set `S` if `BF01` exceeds a chosen threshold `α` (equivalent to a minimum posterior edge probability).
4. Once no adjacency sets are large enough to test further, orient v-structures and propagate orientations with Meek's rules to recover the CPDAG.

## Results

### Simulation studies

The Bayesian and frequentist approaches were compared over:
- Number of nodes `q ∈ {50, 200, 1000}`
- Sample size `n ∈ {100, 200, 500}`
- Sparse random DAGs with edge probability `w = 2/q`
- Frequentist significance level `α = 0.01`; Bayesian posterior threshold `p_th = 0.9`

Evaluated via Structural Hamming Distance (SHD), precision, and recall. The Bayesian method consistently achieved **lower SHD** than the frequentist PC algorithm, with the gap widening substantially as `q` grows to 1000, where the frequentist approach becomes structurally unstable.

### Real data: Riboflavin dataset

Applied to the riboflavin (vitamin B2) production dataset (gene expression in *Bacillus subtilis*), screened down from 4088 to the top 1000 genes by empirical variance (`q = 1000 ≫ n = 71`).

- The Bayesian method recovered a richer skeleton (1053 edges vs. 513 for the frequentist method), while still containing 97% of the frequentist edges — evidence that the frequentist test loses power at this sample size rather than being more conservative.
- Local analysis around the `RIBO` node identified two direct predictors, `YDDK` and `YCKE`, with strong posterior edge probability.
- A global hub analysis identified `YOPZ` as the most connected gene (max degree 5), consistent with a sparse, decentralized regulatory topology rather than a small number of global master regulators.

## Repository structure

> Update this section to match the actual repository layout.

```
.
├── R/ or src/         # Implementation of the Bayesian PC algorithm
├── simulations/        # Simulation study code and results
├── data/                # Riboflavin dataset / preprocessing scripts
├── presentation/        # Slides (BeraldoGattaMarazNaccarelliPoliPruenster.pdf)
└── README.md
```

## References

- Bühlmann, P., Kalisch, M., & Meier, L. (2014). High-Dimensional Statistics with a View Toward Applications in Biology. *Annual Review of Statistics and Its Application*, 1, 255–278.
- Castelletti, F., & Consonni, G. (2024). Bayesian Sample Size Determination for Causal Discovery. *Statistical Science*, 39(2), 305–321.
- Castelletti, F., & Mascaro, A. (2022). BCDAG: An R package for Bayesian structure and Causal learning of Gaussian DAGs. *arXiv:2201.12003*.
- Kalisch, M., & Bühlmann, P. (2007). Estimating High-Dimensional Directed Acyclic Graphs with the PC-Algorithm. *Journal of Machine Learning Research*, 8, 613–636.
- Kalisch, M., Mächler, M., Colombo, D., Maathuis, M. H., & Bühlmann, P. (2012). Causal Inference Using Graphical Models with the R package pcalg. *Journal of Statistical Software*, 47(11), 1–26.
- Spirtes, P., Glymour, C., & Scheines, R. (2000). *Causation, Prediction, and Search* (2nd ed.). MIT Press.
- Tsagris, M. (2019). Bayesian Network Learning with the PC Algorithm: An Improved and Correct Variation. *Applied Artificial Intelligence*, 33(2), 101–123.


## 💻 Repository Workflow Guide

This guide outlines the standard Git workflow for collaborating on this repository, focusing on using a Personal Access Token (PAT) for authentication (required by GitHub over HTTPS) and the essential three-step cycle: Stage, Commit, Push.

### 1. Authentication Setup (Using a Personal Access Token)

- GitHub requires a Personal Access Token (PAT) instead of your regular password for all command-line operations (like git push) over HTTPS.

#### 🔑 A. Generate the PAT

- Go to GitHub: Navigate to Settings > Developer settings > Personal access tokens > Tokens (classic).

- Generate a New Token: Click Generate new token (classic).

- Define scope: You must check the repo box. This grants the token read/write access to your repositories.

- Save the Token: Once generated, IMMEDIATELY COPY the resulting string. This is the only time you will ever see it.

#### 💾 B. Use and Store the PAT

The first time you run a command that requires authentication, Git will prompt you for credentials:

- Username: Enter your GitHub username.

- Password: PASTE THE PAT (the long string) here, NOT your actual GitHub account password.

Note: *Credentials will be stored, but if they are needed, repeat the process*


### 2. Initial Setup: Cloning the Repository

Cloning downloads a copy of the repository to your local computer and sets up the necessary links to GitHub.

- Navigate: Open your terminal or Git Bash and move to the directory where you want to store the project: ```cd ~/Projects/```

- Clone: Copy the HTTPS or SSH link from the GitHub repository page and clone it: ```git clone [repository_url]```

- Enter Directory: ```cd [repository-name]```


### 3. Daily Workflow: Making and Pushing Changes

The core Git workflow is a three-step cycle that turns local file modifications into recorded history on GitHub.

#### Step 1: Check Status and Modify Files

- Download the latest changes from GitHub (*useful to use before you start*): ```git pull``` 

- Shows which files are modified (edited), untracked (new), or staged (ready to commit): ```git status```

#### Step 2: Stage and Commit (Saving the Change)

You must tell Git which changes you want to include in the next version snapshot.

##### 📁 A. Adding / Editing a File
To include new files or modified files in your commit, you must Stage them first:

- To stage a single file: ```git add my_new_file.txt```

- To stage ALL modified, new, and deleted files: ```git add```

##### ✍️ B. Committing the Snapshot
Once files are staged, you create a permanent snapshot (a commit) in your local history:

```git commit -m "Brief but descriptive summary of the changes"```

##### 🗑️ C. Deleting a File
If you want to remove a file from the repository, you must stage the deletion:

- Deletes the file locally AND stages the deletion for the next commit: ```git rm file_to_delete.txt ```

- To commit the deletion: ```git commit -m "Removed old file_to_delete.txt"```

#### Step 3: Push to Remote (Updating GitHub)
After working on the repository, you need to push (upload) all the new commits from your local directory to GitHub.

- To push to Github repository: ```git push```

Note: *If this is the first time you've pushed a new branch, you may need to use:* ```git push -u origin [branch-name]```

Note: *If your file is present locally but not on GitHub, you likely missed the git add step. Always run git status to verify your file is in the "Changes to be committed" (green) section before committing.*
