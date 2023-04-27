#
# experiment.R
#
# run simulation experiments
#

source('../helper_function.R')
source('../mom_ls.R')
source('setting_helper.R')
source('data_helper.R')
source('tuning_helper.R')
source('evaluation_helper.R')
source('bMIND_helper.R')
library(dplyr)
library(nnls)

# -
# 1. set parameters for the experiment
# -

# take arguments from command line
library("optparse")
option_list = list(
  make_option(c("--n_rep"), type="integer", default=1,
              help="number of simulation replications", metavar="integer"),
  make_option(c("--n"), type="integer", default=150,
              help="number of bulk samples", metavar="integer"),
  make_option(c("--p"), type="integer", default=100,
              help="number of genes", metavar="integer"),
  make_option(c("--K"), type="integer", default=2,
              help="number of cell types", metavar="integer"),
  make_option(c("--log_var"), type="numeric", default=8.0,
              help="log variance", metavar="numeric"),
  make_option(c("--cor_model"), type="character", default="MA",
              help="AR, MA or dcSBM for correlated gene clusters", metavar="character"),
  make_option(c("--n_val"), type="integer", default=150,
              help="size of the independent validation samples", metavar="numeric"),
  # parameters when K=2
  make_option(c("--beta1"), type="integer", default=1,
              help="1st parameter for the Beta distribution", metavar="numeric"),
  make_option(c("--beta2"), type="integer", default=2,
              help="2nd parameter for the Beta distribution", metavar="numeric"),
  # parameters specific to AR(1) and MA(1) simulations
  make_option(c("--rho1"), type="numeric", default=0.9,
              help="AR parameter for cell type 1", metavar="numeric"),
  make_option(c("--rho2"), type="numeric", default=0.9,
              help="AR parameter for cell type 2", metavar="numeric"),
  # parameters specific to K=4 experiments
  make_option(c("--equal_strength"), type="character", default="T",
              help="whether cell type 4 have equal expression levels as 1-3", metavar="character"),
  # parameters specific to sensitivity analysis
  make_option(c("--sensitivity"), type="character", default='F',
              help='T or F: whether to add noise to pi', metavar="character"),
  make_option(c("--kappa"), type="numeric",
              default=0, help="scale for the pi noise",
              metavar='numeric'),
  make_option(c("--b"), type="numeric",
              default=0, help="relative magnitude of the bias",
              metavar='numeric'),
  # parameters specific to which LS methods to use
  make_option(c("--method_var"), type="character", default="ols",
          help="run ols or nnls for estimating cell-type-specificvariances", 
          metavar="character"),
  make_option(c("--method_covar"), type="character", default="ols",
          help="run ols or irls for estimating cell-type-specificvariances", 
          metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## parameters specific to simulation experiments
# number of replications
n_rep <- opt$n_rep
# sample seeds for random samples
seed <- 1
set.seed(seed)
seeds <- sample.int(100000, n_rep)

## parameters for the simulation setting
# number of genes in the network
p <- opt$p
# sample size
n <- opt$n
# number of cell types
K <- opt$K
# log variance
log_var <- opt$log_var
# correlation model in ct-specific cluster
cor_model <- opt$cor_model
if(cor_model %in% c('AR', 'MA')){
	# AR(1) strength
	rhos <- c(opt$rho1, opt$rho2)
}
if(K == 2){
	# parameters for beta distribution, s.t. beta ~ beta(beta_1, beta_2)
	beta1 <- opt$beta1
	beta2 <- opt$beta2
	betas <- c(beta1, beta2)
}

## parameters for CSNet regressions
methods <- list(var = opt$method_var,
	covar = opt$method_covar)

## parameters for sensitivity analysis
# whether to conduct sensitivity analysis
sensitivity <- opt$sensitivity == 'T'
if(sensitivity){
  kappa <- opt$kappa
  b <- opt$b
}

## parameters specific to K=4 analysis
# whether to specify equal strength
equal_strength <- opt$equal_strength == 'T'

## fixed parameters

# set S=1 million s.t. simulated counts can be interpreted as FPKM
# (not considering the heterogeneity in gene lengths)
# set Tsize to be a smaller integer
# to decrease the level of bias in simulated correlations
NB_exper_pars <- list(Tsize=200, S=6e+07)

cv_mode <- 'independent'
cor_struct <- 'three_clusters'
gen_th_op <- 'SCAD'
n_val <- 150 # number of samples for independent validation
exp_model <- 'NB'

print(log_var)
print(equal_strength)
# set file prefix
if(K == 2){
	if(cor_model %in% c('AR', 'MA')){
		cor_model_prefix <- sprintf('%s_rho_%.1f_%.1f', cor_model, rhos[1], rhos[2])
	}else{
		cor_model_prefix <- 'dcSBM'
	}
	occ_prefix <- sprintf('_log_var_%.1f_equal_strength_%s', log_var, equal_strength)
	prefix <- sprintf('%s/K_%i/var_%s_covar_%s%s/n_%i_p_%i', 
		cor_model_prefix, K, methods$var, methods$covar, 
		ifelse((log_var != 8.0) | (!equal_strength), occ_prefix, ''),
		n, p)
}

if(sensitivity){
	prefix <- sprintf('%s/%s/kappa_%.2f_b_%.2f', sensitivity, prefix, kappa, b)
}

fig_prefix <- sprintf('figures/%s', prefix)
result_prefix <- sprintf('results/%s', prefix)

print(fig_prefix)
print(result_prefix)

if(!file.exists(file.path(result_prefix))){
  dir.create(fig_prefix, recursive = T)
  dir.create(result_prefix, recursive = T)
}

# -
# 2. simulate gene expression data according to specific parameters 
# -

# generate the simulation setting
sim_setting <- gen_sim_setting(cor_model,
	n, p, log_var, 
	K, betas, 
	AR_MA_ctrl = list(rhos = rhos),
	K_4_ctrl = list(equal_strength = equal_strength),
	NB_exper_pars = NB_exper_pars, 
	verbose=F)

# simulate expression data & run methods
# keep track of errors
coexp_methods <- c('Bulk', 'd-CSNet', 'CSNet', 'bMIND', 's-bMIND', 'bMIND-noninf', 's-bMIND-noninf', 'ENIGMA', 's-ENIGMA')
metrics <- c('F_norm', 'op_norm', 'FPR', 'TPR')
error_mat <- matrix(NA, nrow = n_rep, ncol = K)
error_array <- array(NA, dim = c(n_rep, 4, K), 
  dimnames = list(paste0('rep_', 1:n_rep),
    metrics, paste0('cell_type_', 1:K)))

error_rec <- lapply(1:length(coexp_methods), function(i) error_array)
names(error_rec) <- coexp_methods

for(i_rep in 1:n_rep){
  # simulate cell-type-specific expression
  train_list <- sim_exp(cor_model, n, seeds[i_rep], sim_setting, props = NULL, verbose=F)
  valid_list <- sim_exp(cor_model, n, seeds[i_rep]*2 + 10, sim_setting, 
    props = sim_setting$props[1:n_val, ],verbose=F)
  
  # evaluate ct-specific co-expressions estimated by different methods
  

  # -
  # evaluate CSNet
  # -

  # obtain d-CSNet estimates
  train_cov_est <- mom_ls(train_list$data$P, train_list$data$X, methods)
  train_cor_est <- lapply(train_cov_est, function(x) get_cor_from_cov(x)) 
  error_rec[['d-CSNet']][i_rep, , ] <- sapply(1:K, function(k){
    eval_errors(train_cor_est[[k]], sim_setting$R[[k]], metrics, dense_estimate = T)
  })

  # tune thresholding parameter with the independent validation data
  tuning_result <- th_tuning(train_list$data, valid_list$data,
    coexp_method = 'CSNet',
    ctrl_list = list(ls_methods = methods))

  # obtain CSNet estiamtes
  train_cor_th_est <- lapply(1:K, function(k){
    generalized_th(train_cor_est[[k]], tuning_result$th[k], gen_th_op, F)
  })
  error_rec[['CSNet']][i_rep, , ] <- sapply(1:K, function(k){
    eval_errors(train_cor_th_est[[k]], sim_setting$R[[k]], metrics)
  })

  # -
  # evaluate sparse bulk estimate
  # -
  bulk_train_cor_est <- cor(train_list$data$X)
  bulk_tuning_result <- th_tuning(train_list$data, valid_list$data, 
    coexp_method = 'Bulk')
  bulk_train_cor_th_est <- lapply(1:K, function(k){
    generalized_th(bulk_train_cor_est, bulk_tuning_result$th, gen_th_op, F)
  }) # same for two cell types
  error_rec[['Bulk']][i_rep, ,] <- sapply(1:K, function(k){
    eval_errors(bulk_train_cor_th_est[[k]], sim_setting$R[[k]], metrics)
  })

  # -
  # evaluate bMIND
  # -
  
  ## with informative prior
  # bMIND
  bMIND_train_est <- coexp_by_bMIND(train_list, prior_info = 'informative')
  error_rec[['bMIND']][i_rep, , ] <- sapply(1:K, function(k){
    eval_errors(bMIND_train_est[[k]], sim_setting$R[[k]], metrics)
  })

  # sparse bMIND
  bMIND_tuning_result <- th_tuning(train_list, valid_list,
    coexp_method = 'bMIND',
    ctrl_list = list(prior = 'informative'))
  sparse_bMIND_train_th_est <- lapply(1:K, function(k){
    generalized_th(bMIND_train_est[[k]], bMIND_tuning_result$th[k], gen_th_op, F)
  })
  error_rec[['s-bMIND']][i_rep, , ] <- sapply(1:K, function(k){
    eval_errors(sparse_bMIND_train_th_est[[k]], sim_setting$R[[k]], metrics)
  })

  ## with non-informative prior
  # bMIND
  bMIND_noninf_train_est <- coexp_by_bMIND(train_list, prior_info = 'non-informative')
  error_rec[['bMIND-noninf']][i_rep, , ] <- sapply(1:K, function(k){
    eval_errors(bMIND_noninf_train_est[[k]], sim_setting$R[[k]], metrics)
  })

  # sparse bMIND
  bMIND_noninf_tuning_result <- th_tuning(train_list, valid_list,
    coexp_method = 'bMIND',
    ctrl_list = list(prior = 'non-informative'))
  sparse_bMIND_noninf_train_th_est <- lapply(1:K, function(k){
    generalized_th(bMIND_noninf_train_est[[k]], bMIND_noninf_tuning_result$th[k], gen_th_op, F)
  })
  error_rec[['s-bMIND-noninf']][i_rep, , ] <- sapply(1:K, function(k){
    eval_errors(sparse_bMIND_noninf_train_th_est[[k]], sim_setting$R[[k]], metrics)
  })

  # -
  # evaluate ENIGMA
  # -

  # est <- mom_ls(sim_setting$props, bulk_exp, 
  # methods = methods)
  # saveRDS(est, 'new_est.rds')
}
# saveRDS(sim_setting, 'new_setting.rds')
# saveRDS(gene_exp_list, 'new.rds')

for(coexp_m in coexp_methods[1:3]){
  cat('\n')
  print(coexp_m)
  for(m in metrics){
    print(m)
    for(k in 1:K) print(sprintf('cell type %i: %.2f (%.2f)', k, mean(error_rec[[coexp_m]][, m, k]), sd(error_rec[[coexp_m]][, m, k])))
  }
}



