#!/bin/bash

# Steps to run CIBERSORTx
# 1. Contrust single cell reference matrix
# 2. Set up CIBERSORTx
# 3. Run CIBERSORTx

# 1. Construct single cell reference matrix
Rscript make_sc_input.R

# 2. Set up CIBERSORTx
#   2.1 Apply for an account and a token at https://cibersortx.stanford.edu 
#   2.2 Download the Docker implementation at https://hub.docker.com/r/cibersortx/hires

# 3. Run CIBERSORTx with singularity
# username and token were generated in step 2
singularity exec \
--bind ${data_dir}:/src/data,${data_dir}:/src/outdir \
cbx_highres.sif \
/src/CIBERSORTxFractions \
--username ${username} --token ${token} \ 
--single_cell TRUE \
--refsample output/ROSMAP_sc_ref_all_group_5000.txt \
--mixture bulk_train.txt \
--fraction 0 \
--rmbatchSmode TRUE \
--outdir ${outdir} 
