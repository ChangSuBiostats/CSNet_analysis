library(dplyr)
library(optparse)
source('table_helper.R')

option_list = list(
  make_option(c("--result_prefix"), type="character", 
              default="results/MA_rho_0.5_0.5/K_2/log_var_8.0_equal_strength_TRUE/n_150_p_100",
              help="which experiment to print", metavar="character"),
  make_option(c("--setting"), type="character",
              default="standard",
              help="which set of methods to print", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
result_prefix <- opt$result_prefix
setting <- opt$setting

i_ind <- ifelse(grepl('sensitivity', result_prefix), 1, 0)
n_p_setting <- strsplit(result_prefix, '/')[[1]][5 + i_ind]
n <- strsplit(n_p_setting, '_')[[1]][2] %>% as.integer
p <- strsplit(n_p_setting, '_')[[1]][4] %>% as.integer
K <- strsplit(strsplit(result_prefix, '/')[[1]][3 + i_ind], '_')[[1]][2] %>% as.integer
if(grepl('sensitivity', result_prefix)){
  kappa_b_setting <- strsplit(result_prefix, '/')[[1]][7]
  kappa <- strsplit(kappa_b_setting, '_')[[1]][2] %>% as.numeric
  b <- strsplit(kappa_b_setting, '_')[[1]][4] %>% as.numeric
}

n_rep <- 200
error_rec_1 <- readRDS(sprintf('%s/n_rep_%i_i_rep_%i_error.rds', result_prefix, n_rep, 1))
error_mat <- matrix(NA, nrow = n_rep, ncol = nrow(error_rec_1[[1]]))
error_list <- lapply(1:length(error_rec_1), function(i) lapply(1:K, function(k) error_mat))
names(error_list) <- names(error_rec_1)

for(i_rep in 1:n_rep){
  tmp <- readRDS(sprintf('%s/n_rep_%i_i_rep_%i_error.rds', result_prefix, n_rep, i_rep))
  for(m in names(error_list)){
    for(k in 1:K){
      error_list[[m]][[k]][i_rep,] <- tmp[[m]][,k]
    }
  }
}

error_methods_names <- names(error_list)
if('d-oracle' %in% error_methods_names){
  error_methods_names['oracle' == error_methods_names] <- 's-oracle'
  error_methods_names['d-oracle' == error_methods_names] <- 'oracle' 
  names(error_list) <- error_methods_names
}


if(setting == 'standard'){
  coexp_methods <- c('Bulk', 'd-CSNet', 'CSNet', 'bMIND', 's-bMIND')
}else if(setting == 'all'){
  coexp_methods <- c('Bulk', 'd-CSNet', 'CSNet', 'd-CSNet-ols', 'CSNet-ols', 'bMIND', 's-bMIND', 'bMIND-noninf', 's-bMIND-noninf', 'ENIGMA', 's-ENIGMA')#, 'ENIGMA-trace', 's-ENIGMA-trace')
}else if(setting == 'reproduce'){
  coexp_methods <- c('Bulk', 'd-CSNet-ols', 'CSNet-ols', 'bMIND')
}else if(setting == 'ols_vs_wls'){
  coexp_methods <- c('d-CSNet-ols', 'd-CSNet', 'CSNet-ols', 'CSNet')
}else if(setting == 'bMIND-inf'){
  coexp_methods <- c('bMIND-inf','s-bMIND-inf')
}else if(setting == 'ENIGMA'){
  coexp_methods <- c('d-CSNet', 'ENIGMA', 'CSNet', 's-ENIGMA')
}else if(setting == 'oracle'){
  coexp_methods <- c('d-CSNet', 'oracle', 'CSNet', 's-oracle')
}else if(setting == 'sparse_bMIND'){
  coexp_methods <- c('d-CSNet', 'bMIND', 'CSNet', 's-bMIND')
}else if(setting == 'sensitivity_CSNet'){
  coexp_methods <- c('CSNet')
}else if(setting == 'sensitivity_bMIND'){
  coexp_methods <- c('s-bMIND')
}else if(setting == 'small_p'){
  coexp_methods <- c('d-CSNet', 'bMIND', 'ENIGMA', 'CSNet', 's-bMIND', 's-ENIGMA')
}

if(grepl('sensitivity', result_prefix)){
  # evaluate the difference between noisy proportions and true proportions
  sim_setting <- readRDS(sprintf('%s/sim_setting.rds', result_prefix))
  noisy_pi <- sim_setting$props
  truth_prefix <- 'results/AR_10_rho_0.8_0.8/K_2/log_var_8.0_equal_strength_TRUE/n_600_p_100'
  true_pi <- readRDS(sprintf('%s/sim_setting.rds', truth_prefix))$props
  
  rho_cor <- cor(noisy_pi[,1], true_pi[,1])
  rMSE <- sqrt(mean((noisy_pi - true_pi)^2)) / mean(true_pi)
  
  print_a_sensitivity_setting(error_list, coexp_methods, kappa, b, rho_cor, rMSE)
}else{
  print_a_setting(error_list, coexp_methods, n, p)
}
#cat('\\\\\\hline')
