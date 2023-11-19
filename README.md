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
The analysis was conducted in R and Conda was used as the package environment tool. 'requirement.txt' specifies the Conda environment for reproducing the analysis in this paper. We note that this list is very long, and a shorter list of key required packages are:

* Network estimation: nnls (1.4), MIND (0.3.3, installed from GitHub [repo](https://github.com/randel/MIND)), ENIGMA (0.1.6, installed from GitHub [repo](https://github.com/WWXkenmo/ENIGMA))
* Tools:
  * Data manipulation: dplyr (1.0.10), reshape2 (1.4.4), Matrix (1.5-4)
  * Software installation from GitHub: devtools (2.4.5)
  * Taking parameters for R scripts: optparse (1.7.3)
  * Visualization: WGCNA (1.71), ggplot2 (3.4.2), gridExtra (2.3), grid (4.2.2), ggpubr (0.6.0), pheatmap (1.0.12), RColorBrewer (1.1-3), ROCR (1.0-11), scales (1.2.1)
  * Sensitivity analyses: DirichletReg (0.7-1)
* Real data analysis: 
  * Deconvolution: CIBERSORTx implementation on Docker [link](https://hub.docker.com/r/cibersortx/hires)
  * Annotating gene names based on ENSEMBL ID: biomaRt (2.54.1)

# data
This folder saves the real data correlation structure used in simulation, and the list of four gene sets used in real data analysis. See 'data/README' for more details.

