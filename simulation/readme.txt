Simulation Setups

We generated variable timetrees from a birth-death process using TreeSim in R (Stadler, 2011), with a birth rate of 5.0 and a death rate of 4.0, conditioned on a root age of 1.0. The ages are on a relative scale and can be arbitrarily rescaled depending on the chosen time unit. We kept 100 trees with no more than 50 tips to make sure the data size is manageable.

For each tree, we simulated the evolution of continuous traits under a BM process with a diffusion rate of 1.0. The starting state at the root was set to 0.0. We recorded a moderate size of 200 continuous traits in the data matrix. In addition to the complete data, we also generated datasets with 50% missing states in the extinct taxa and 10% missing in the extant taxa, mimicking observations in empirical datasets. Specifically, we replaced each state by a question mark with the corresponding probability.

For comparisons, we also generated 200 variable binary discrete characters under the Mkv model (Lewis, 2001). The substitution rate was set to 1.0. And to test the functionality of supporting mixed data types, we generated data matrices each containing 100 continuous and 100 discrete characters.

We performed both nonclock and tip-dating analyses for each dataset using MrBayes (Ronquist et al., 2012). In the nonclock analyses, we used uniform prior for the tree topologies and gamma-Dirichlet prior for the branch lengths (Zhang et al., 2012). In the tip-dating analyses, we used the constant-rate fossilized birth-death (FBD) prior (Stadler, 2010) for the timetrees and independent-lognormal relaxed clock model (Drummond et al., 2006) for the evolutionary rates (although the true model is strict clock). For each inference, two independent chains were executed for enough generations to ensure convergence and sufficient effective sample sizes (ESS > 100). The posterior tree samples are summarized as a 50% majority-rule consensus tree.

Citations

Stadler, T. (2011). Simulating trees with a fixed number of extant species. Systematic Biology, 60(5), 676–684.
Lewis, P. O. (2001). A likelihood approach to estimating phylogeny from discrete morphological character data. Systematic Biology, 50(6), 913–925.
Ronquist, F., Teslenko, M., Mark, P. van der, Ayres, D. L., Darling, A., Höhna, S., Larget, B., Liu, L., Suchard, M. A., & Huelsenbeck, J. P. (2012). MrBayes 3.2: efficient Bayesian phylogenetic inference and model choice across a large model space. Systematic Biology, 61(3), 539–542.
Zhang, C., Rannala, B., & Yang, Z. (2012). Robustness of compound Dirichlet priors for Bayesian inference of branch lengths. Systematic Biology, 61(5), 779–784.
Stadler, T. (2010). Sampling-through-time in birth-death trees. Journal of Theoretical Biology, 267(3), 396–404.
Drummond, A., Ho, S., Phillips, M., & Rambaut, A. (2006). Relaxed Phylogenetics and Dating with Confidence. PLoS Biology, 4(5), e88.
