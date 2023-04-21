# A joint Bayesian hierarchical model for estimating SARS-CoV-2 diagnostic and subgenomic RNA viral dynamics and seroconversion

Here are some example R code and JAGS code for the model implementation in "A joint Bayesian hierarchical model for estimating SARS-CoV-2 diagnostic and subgenomic RNA viral dynamics and seroconversion" by Dong and Brown (2022). 

1. *joint_vl_sero_r.R*: R code for fitting the joint Bayesian model. 

2. *joint_vl_sero_jags.txt*: JAGS code for fitting the joint Bayesian model. 

3. *model_par_diag_primary.pdf*: Full model diagnostic report for the primary model fitted to the COVID-19 PEP study data, as described in Section 4.1 of the main manuscript and Section S1.3 of the Supplement. The *ggmcmc* package in R (version 1.5.1.1) was used to create this report. Specifically, we examined the following plots and summary metrics for key model parameters: 

- Histograms of posterior samples.
- Density plots with different colors by chain. 
- Traceplots with different colors by chain to assess convergence and chain problems.
- Plots of the running means to check how quickly the chain is approaching its target distribution. 
- Overlapped density plots that compare the last 10 percent of the chain with the whole chain. 
- Autocorrelation plots to check for autocorrelation of posterior samples. 
- Crosscorrelation plot to diagnose potential problems of convergence due to highly correlated parameters.
- The Potential Scale Reduction Factor ($\hat{R}$)(\citealp{Gelman2013BayesianDA}) for comparison of the between-chain variation with the within-chain variation. 
- The Geweke z-score diagnostic (\citealp{geweke1992evaluating}) for comparison of the first part of the chain with its last part. 
