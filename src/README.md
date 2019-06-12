# src folder

This folder contains the code we wrote to process and analyze our data. The analysis starts from the raw collections data from 2015, found in `raw_data/erato_collected_2015.csv`.

0. Throughout our analysis, we wrote and re-used a number of functions that are used in multiple other scripts. They are found in the `functions.R` and `hzar_functions.R` files. 

1. First, we checked that the subsites at which we sampled did not show significant differentiation. That was done in `subsite_differentiation.R`, and revealed that all subsites could be merged. 

2. Then, we combined sampling data from all three timepoints into a single transect using the `joint_transect.R` script. 

3. Next, we calculated the allele frequencies for the three alleles segregating in some populations, using `three_allele_freqs.R`. 

4.1 We then developed and tested new cline models. The `Stan` code for those models is found in the `stan_models` folder, and there is a separate file for each of the five possible introgression tail models.

4.2 We simulated 720 datasets and then analyzed each using both our Bayesian method and a maximum likelihood model in the R package `HZAR`. The code to run those simulations is in `simulate_bayes_HZAR.R`, while the code to compare the results of the two models is in `compare_bayes_HZAR.R`

5. We fit clines using our Bayesian approach for each of the three timepoints, using the `cline_multi_mallet.R`, `cline_multi_blum.R`, and `cline_multi_thurman.R` scripts.  

6. We performed model checks on the cline models using the `cline_multi_modelchecks.R` script, and then analyzed the results with the `cline_multi_resuls.R` script. 

7. We did environmental analysis of forest cover using the `forest_cover.R` script. 

8. Finally, we generated figures for publication using the `figure_1.R`, `figure_2.R`, and `figure_3.R` scripts. 



