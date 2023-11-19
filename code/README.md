This folder contains all codes for reproducing the numerical results in the manuscript.

* mom_ls.R and cross_validation.R are the implementation of the methodology described in Section 2 of the paper.

* The folder simulation/ contains the workflow (README.md) and codes for simulation experiments in Section 3.

* The folder real_data/ contains the workflow (README.md) and codes for real data experiments in Section 4.

* 'requirement.txt' specifies the Conda environment for reproducing the analysis in this paper. We note that this list is very long, and a shorter list of key required packages are:

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
