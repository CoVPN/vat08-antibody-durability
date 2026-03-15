## Summary of sharing of data and code for Li et al., "Durability of Neutralizing and Anti-Spike Binding IgG Antibody Responses to Monovalent (D614) and Bivalent (D614 + B.1.351) AS03-Adjuvanted SARS-CoV-2 Recombinant Protein Vaccines."

**Li Li**

**March 14, 2026**

- All analyses are based on the analysis-ready data files `vat08_combined_data_processed_longitudinal_bAb.csv` and `vat08_combined_data_processed_longitudinal_nAb.csv`.
- A single statistical report includes all the main figures and tables. A copy of the statistical report is on the SCHARP network drive `T:\covpn\p3005\analysis\antibody_durability`


### Figure 1
nAb-ID50 titers against Omicron BA.4/BA.5 and bAb-IgG Spike concentrations against Omicron before detected SARS-CoV-2 infection in the nAb and bAb longitudinal analysis cohorts, 
separately by stage, treatment arm, and Naïve vs. Non-naïve status. 

**Statistical report**: vat08_antibody_durability_ms.pdf

- Data were extracted from page 2.

- To reproduce this report, follow the following instructions:

```bash

# a) download the zip file that contains the code from https://github.com/CoVPN/vat08-antibody-durability

# b) install R packages 
# Assume that R 4.4.3 and R package `remotes` have been installed

R
    remotes::install_version("tidyverse", version = "2.0.0")
    remotes::install_version("plyr", version = "1.8.9")
    remotes::install_version("lme4", version = "1.1-37")
    remotes::install_version("doParallel", version = "1.0.17")
    remotes::install_version("doRNG", version = "1.8.6.1")
    remotes::install_version("ggplot2", version = "3.5.2")
    remotes::install_version("scales", version = "1.4.0")

# c) set the directory
R
   setwd("~/vat08-antibody-durability-master")
   repoDir <- here::here()
   
# d) source R files
R 
    source(file.path(repoDir, "code/desc/descFig.R"))
    source(file.path(repoDir, "code/LMMmodeling/confidenceInt.R"))
    source(file.path(repoDir, "code/LMMmodeling/permutation_test.R"))
    source(file.path(repoDir, "code/LMMmodeling/FittedTrajectoryPlot.R"))
    source(file.path(repoDir, "code/summary/summary.R"))
    
# e) generate report pdf

R
    rmarkdown::render(file.path(repoDir, "code/vat08_antibody_durability_ms.Rmd"), 
    output_file = file.path(repoDir, "code/vat08_antibody_durability_ms.pdf"), quiet=TRUE)

```

### Figure 2
Empirical geometric mean (GM) nAb-ID50 titers against BA.4/5 and against Reference and fitted linear mixed models by treatment arm and stage for the Non-naïve groups in the nAb longitudinal analysis cohort.
Estimated D43 GM nAb-ID50 titers, durability, and D202-to-D43 geometric mean ratios (GMR) by treatment arm, with estimated vaccine-to-placebo ratios (95% CI).

**Statistical report**: vat08_antibody_durability_ms.pdf


- Data were extracted from page 3 and 4.

### Figure 3
Empirical geometric mean (GM) bAb-IgG Spike concentrations against Omicron and against the Index strain and fitted linear mixed models by treatment arm and stage for the Non-naïve groups in the bAb longitudinal analysis cohort.
Estimated D43 GM bAb-IgG Spike concentrations, durability, and D202-to-D43 geometric mean ratios (GMR) by treatment arm, with estimated vaccine-to-placebo ratios (95% CI).

**Statistical report**: vat08_antibody_durability_ms.pdf

- Data were extracted from page 5 and 6.

### Figure 4
Empirical geometric mean (GM) nAb-ID50 titers among Non-naive participants and Naive participants and fitted linear mixed models by Naive/Non-naive status and stage in the nAb longitudinal analysis cohort.
Estimated D43 GM nAb-ID50 titers, durability, and D202-to-D43 geometric mean ratios (GMR), with estimated BV-to-MV ratios (95% CI).

**Statistical report**: vat08_antibody_durability_ms.pdf

- Data were extracted from page 7 and 8.

### Figure 5
Empirical geometric mean (GM) bAb-IgG Spike concentrations among Non-naive participants and Naive participants and fitted linear mixed models by Naive/Non-naive status and stage in the bAb longitudinal analysis cohort.
Estimated D43 GM bAb-IgG Spike concentrations, durability, and D202-to-D43 geometric mean ratios (GMR), with estimated BV-to-MV ratios (95% CI).

**Statistical report**: vat08_antibody_durability_ms.pdf

- Data were extracted from page 9 and 10.


### Table 1
Comparisons of estimated geometric mean (GM) D43 antibody marker levels, durability, and D202-to-D43 geometric mean ratio (GMR)s between the Non-naïve vaccine and Naïve vaccine recipients in the nAb and bAb longitudinal analysis cohorts.

**Statistical report**: vat08_antibody_durability_ms.pdf

- Data were extracted from page 11.


