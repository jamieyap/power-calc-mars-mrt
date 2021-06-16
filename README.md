<!---
---
output:
  pdf_document: default
urlcolor: magenta
---
--->

# \Huge Power Calculation for the Mobile-Assistance for Regulating Smoking Micro-Randomized Trial

\Large _Jamie Yap_ 


\large {email: jamieyap@umich.edu; website: https://github.com/jamieyap}


## About this Repository

The material in this repository is a supplement to the manuscript titled `The Mobile-Assistance for Regulating Smoking (MARS) Micro-Randomized Trial Design Protocol' (Nahum-Shani, et al., 2021), submitted for consideration to the Journal of Contemporary Clinical Trials. This repository contains code and documentation for the calculation of power for the Primary Aim and Secondary Aim of the MARS Micro-Randomized Trial.

The work in this repository was funded by NIH grant U01CA229437 and developed at the d3lab: [https://d3lab.isr.umich.edu/](https://d3lab.isr.umich.edu/)

_Keywords: Causal Inference, Intensive Longitudinal Data, Adaptive Interventions, Mobile Health (mHealth), Micro-Randomized Trial, Power Calculation_


## Start Here

| <img height=0 width=350> File <img height=0 width=350> | <img height=0 width=800> Brief Description <img height=0 width=800> |
|:------------------------------------------:|:--------------------------------------------------------------------------------------------------|
| documentation.pdf | Provides more detail about the set-up and displays the results of the power calculation for the Primary Aim and Secondary Aim of the MARS Micro-Randomized Trial.|

## Replicating Material in this Repository

### Obtaining Power Calculation Results

1. Install the required software to perform computations: R version 4.0.5 and the packages bindata version 0.9-20, geeM version 0.10.1, dplyr version 1.0.6

2. Run run-all-calculations.R

   This will produce the following csv files which are collated into a table displaying results for _Simulation Study 1_ (see the folder sim_study_01)

    * dat_all_results_0.8.csv
    * dat_all_results_0.85.csv
    * dat_all_results_0.9.csv
    * dat_all_results_0.95.csv
    * dat_all_results_1.csv

   This will produce the following csv files which are collated into a table displaying results for _Simulation Study 2_ (see the folder sim_study_02)

    * dat_all_results_0.8_0.1.csv
    * dat_all_results_0.8_0.2.csv

### Displaying Power Calculation Results in a PDF

We note that steps 3 and 4 are _optional_ steps

3. Additionally install the following software: the packages rmarkdown version 2.8, knitr version 1.33, kableExtra version 1.3.4

4. Knit the markdown file documentation.Rmd to display documentation.pdf


## References


1. Friedrich Leisch, Andreas Weingessel and Kurt Hornik (2021). bindata: Generation of Artificial Binary Data. R package version 0.9-20. https://CRAN.R-project.org/package=bindata

2. Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2021). dplyr: A Grammar of Data Manipulation. R package version 1.0.6. https://CRAN.R-project.org/package=dplyr

3. McDaniel, L. S., Henderson, N. C., and Rathouz, P. J. (2013) Fast pure R implementation of GEE: application of the Matrix package The R Journal, 5/1:181--187

4. Nahum-Shani, I., Potter, L., Lam, C., Yap, J., Moreno, A., Stoffel, R., Wu, Z., Dempsey, W., Kumar, S., Murphy, S., Rehg, J., Wetter, D. (2021). The Mobile-Assistance for Regulating Smoking (MARS) Micro-Randomized Trial Design Protocol. Submitted to the Journal of Contemporary Clinical Trials

