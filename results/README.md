# results

This folder has figures, tables, and `.Rdata` files  which contain the results of various analyses.

The `bayes_simulation_res.Rdata` and `corr_simulation_res.Rdata` files contain the results from running our Bayesian method and the corrected-ML liklihood approach in HZAR on 720 simulated data sets. Those files were generated with `src/simulate_bayes_HZAR.R`.

For our analysis, we used two `.Rdata` files containing the model fit results for the cline analysis across each of the three timepoints ("mallet" for the 1986 sampling, "blum" for the 1999-2000 sampling, and "thurman" for the 2015 sampling). One file contains the stanfit object for only the "best" model, that is, the one with the lowest WAIC (e.g., `thurman_bestCline.Rdata`). One file contains the stanfit objects of all 5 cline models (e.g., `thurman_allClines.Rdata`), but those files are over 100MB and are too large to host on GitHub. Instead, they can be found in the Dryad repository for this paper (doi:10.5061/dryad.c157v0v). These files were generated with the `src/cline_multi_mallet.R`, `src/cline_multi_blum.R`, and `src/cline_multi_thurman.R` scripts. 

The `triallelic_freqs.csv` file contains the frequencies of the all three possible alleles at each site, and was generated with the `src/three_allele_freqs.R` script.

Finally, there are the .pdf files for figures 1, 2, and 3. `figure_1_map.pdf` was generated with `src/figure_1.R`. That forms the basis for `figure_1_final.svg`. We used Inkscape to add the butterfly art to that file, and then exported it to the pdf file that appears in the manuscript, `figure_1_final.pdf`.

The `.pdf` files for figure 2 and figure 3 were generated with the `src/figure_2.R` and `src/figure_3.R` scripts, respectively. 

