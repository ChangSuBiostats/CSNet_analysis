CSNet analysis
================

This GitHub repository contains the codes for reproducing the analysis in ''Estimating cell-type-specific gene co-expression networks from bulk gene expression data with an application to Alzheimer’s disease'' by Chang Su, Jingfei Zhang and Hongyu Zhao [[bioRxiv]](https://www.biorxiv.org/content/10.1101/2021.12.21.473558v2).

# code
This folder contains the implementation of CSNet (mom_ls.R and cross_validation.R), and codes for reproducing the analysis in Sections 3 and 4 of the paper (simulation studies and real data studies).

## simulation/
All results in simulation studies (Section 3) can be reproduced following the README under 'code/simulation'. 

## real data/
All results in the real data studies (Section 4) can be reproduced following the README under 'code/real_data'.

## requirements.txt
The analysis was conducted in R and Conda was used as the package environment tool. 'requirement.txt' specifies the Conda environment for reproducing the analysis in this paper. We note that this list is very long, and we provide a shorter list of key required packages in 'code/README.md'.

# data
This folder saves the real data correlation structure used in simulation, and the list of four gene sets used in real data analysis. See 'data/README' for more details.

