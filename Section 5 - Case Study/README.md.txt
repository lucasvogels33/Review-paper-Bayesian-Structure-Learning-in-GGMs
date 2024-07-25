This folder contains the code to reproduce the results of Section 5 ("Case study: Human Gene Expression") of the article: 
"Bayesian Structure Learning in Undirected Gaussian Graphical Models: Literature Review with Empirical Comparison"

This folder contains 5 files:
- edge_convergence.R: this file can run a specified algorithm (BDA or SSO) for a specified initial density (between 0 and 1) on the human gene dataset. It ouputs an .Rdata file used by make_plots.R to create Figure 3
- links_convergence.R: this file can run a specified algorithm (BDA or SSO) on the human gene dataset. It ouputs an .Rdata file used by make_plots.R to create Figure 4.
- networks.R: this file can run a specified algorithm (BDA or SSO) for a specified prior density (between 0 and 1) on the human gene dataset. It ouputs an .Rdata file used by make_plots.R to create Figure 5.
- prior_sensitivity.R: this can run a specified algorithm (BDA or SSO) for a specified prior density (between 0 and 1) on the human gene dataset. It ouputs an .Rdata file used by make_plots.R to create Figure S1 and S2 in the supplementary material.
- make_plots.R: this file downloads the .Rdata files created by edge_convergence.R, links_convergence.R, networks.R and/or prior_sensitivity.R and creates figures 3,4,5 (in the article) and S1 and S2 (in the supplementary material)
