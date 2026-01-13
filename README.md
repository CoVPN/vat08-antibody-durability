## Longitudinal analysis in the VAT00008 Trial
*R* code implementing the data analyses used to generated figures and tables in the manuscript Li et al., 
Durability of Neutralizing and Anti-Spike Binding IgG Antibody Responses to Monovalent (D614) and Bivalent 
(D614 + B.1.351) AS03-Adjuvanted SARS-CoV-2 Recombinant Protein Vaccines.


### 1. System Requirements

  Unless otherwise specified, software was tested on macOS Sequoia 15.6.1 running R version 4.4.3.
  
  Main programs as well as specific R packages and R files required by each program are listed below.

* `code/desc/descFig.R`: generates antibody trajectory plots before and after SARS-CoV-2 infections 
    + tidyverse 2.0.0 
    + plyr 1.8.9 
    + common.R 
    + descFigUtils.R
     
* `code/LMMmodeling/confidenceInt.R`: computes bootstrap confidence intervals for Day 43 Geometric mean, durability, and D202-to-D43 Geometric mean ratio, as 
well as ratios of these metrics for comparisons such as vaccine vs. placebo, Nonnaive vs. naive, Stage 2 vs. Stage 1.
    + tidyverse 2.0.0 
    + plyr 1.8.9 
    + lme4 1.1-37
    + doParallel 1.0.17
    + doRNG 1.8.6.1
    + common.R
    + utils.R
    
    
* `code/LMMmodeling/permutation_test.R`: computes p-values from permutation tests for comparing D43 Geometric means, durability, and D202-to-D43 Geometric mean ratios
    + tidyverse 2.0.0 
    + plyr 1.8.9 
    + lme4 1.1-37
    + doParallel 1.0.17
    + doRNG 1.8.6.1
    + common.R
    + utils.R

    
* `code/LMMmodeling/FittedTrajectoryPlot.R`: generate plots for empirical geometric mean (GM) antibody markers and point-wise GM of fitted responses
    + tidyverse 2.0.0 
    + plyr 1.8.9 
    + lme4 1.1-37
    + ggplot2 3.5.2
    + common.R
    + utils.R

* `code/summary/summary.R`: generate plots and tables that summarize the results of the bootstrap confidence intervals and permutation tests
    + tidyverse 2.0.0 
    + plyr 1.8.9 
    + ggplot2 3.5.2
    + common.R
    + summaryUtils.R

*  `code/MainFiguresTables.Rmd`: compile main tables and figures in a pdf   

### 2. Installation Guide
  
* Install required version of *R*.  
* Install required *R* packages.  
* Clone this repository.
  
### 3. User Instructions
  Input datasets include `vat08_combined_data_processed_longitudinal_bAb.csv`, `vat08_combined_data_processed_longitudinal_nAb.csv`. 
  Place both datasets in the `data` subfolder within the cloned repository. 
  
  Create a new R project for the cloned repository. From the command line, starting in the `code` directory, run the following commands. 
  Each code file should complete in under 10 minutes except for `code/LMMmodeling/confidenceInt.R` and 
  `code/LMMmodeling/permutation_test.R` which may take a few hours, depending on the number of nodes used for parallel computation.
  All output are saved as either figures (`.pdf` files in the `figures` directory) or tables (`.csv` files in the `tables` directory).
    
    R CMD BATCH desc/descFig.R &
    R CMD BATCH LMMmodeling/confidenceInt.R &
    R CMD BATCH LMMmodeling/permutation_test.R &
    R CMD BATCH LMMmodeling/FittedTrajectoryPlot.R &
    R CMD BATCH summary/summary.R 
    R CMD BATCH MainFiguresTables.Rmd
    
    