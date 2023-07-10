#!/bin/bash
#SBATCH --job-name=roc_AR_10
#SBATCH --time=1:00:00

cd ..
ml miniconda
conda activate r_csnet_3

for n in 150 600
do
  for p in 100 200
  do 
    for setting in full #standard
    do
      Rscript plot_ROC.R \
    --result_prefix results/AR_10_rho_0.8_0.8/K_2/log_var_8.0_equal_strength_TRUE/n_"$n"_p_"$p" \
    --setting "$setting"
    done
  done
done
