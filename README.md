# Career Modeling with Missing Data and Traces public repository organization

This describes the organization of scripts within the repository attached to the article "Career Modeling with Missing Data and Traces", Voldoire, Ryder, Lahfa. 

In "data/clean" folder:

- "dataTot.rds" is the data made available for this paper, after anonymization and censoring.

In "scripts/inference" folder:

- "inference_functions.R" is the main script where inference functions are implemented.
- "inference_mcmc.R" is the main function to call to estimate the original models in study 1 and study 2.
- "descriptive_statistics.R" aims at constructing cardinality statistics in both studies.
- "RR_feedback.R" implements all the suggested elements by the reviewers, and is organized by issues.
- 
These scripts construct result data which are only available on Zenodo and should be downloaded before replicating the remaining of the work. 

Final results data scripts available on Zenodo 

- "results/ena_config_X.rds": result data for models X = 1-5  in study 1
- "results/groups_config_X.rds": result data for models X = 1, 2 in study 2
- Within the "revise_resubmit/" folder (see RR_feedback.R for more information):
  
	- "ena_prior.rds": analysis with another prior 
	- "phi_XX.rds": analysis of study 1, substantive model with different weights for phi
	- "phi_all.rds": combination of multiple weights (used to study correlation across phis)
	- "syn_power.rds": study on synthetic data, showing the data augmentation scheme
	- "syn_miss.rds": study on synthetic data, showing that standard HMM behaves poorly when misspecified (motivating Markov switching)
	- "ena_CVX.rds": cross-validation analysis of reconstruction loss for models in study 1.
	- "groups_CVX.rds": cross-validation analysis of reconstruction loss for models in study 2.
	- "convergence_analysis.rds": replication of models in study 1 to compute some Rhat on multiple chains. 


Other scripts:
- "scripts/inference/tables_crafter.R" crafts tables and graphs that are present in the paper (can only be run with downloaded final data from Zotero)
