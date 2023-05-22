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
source('ENIGMA_helper.R')
source('visualization_helper.R')
source('sensitivity_helper.R')
library(dplyr)
library(nnls)
library(ggplot2)
library(gridExtra)

# -
# 1. set parameters for the experiment
# -

# take arguments from command line
library("optparse")
option_list = list(
  make_option(c("--n_rep"), type="integer", default=1,
              help="number of simulation replications", metavar="integer"),
  make_option(c("--i_rep"), type="integer", default=1,
              help="the simulation replication", metavar="integer"),
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
  # make_option(c("--method_var"), type="character", default="ols",
  #         help="run ols or nnls for estimating cell-type-specificvariances", 
  #         metavar="character"),
  # make_option(c("--method_covar"), type="character", default="ols",
  #         help="run ols or irls for estimating cell-type-specificvariances", 
  #         metavar="character"),
  # whether to save estimates
  make_option(c("--save_est"), type="character", default="F",
          help="whether to save co-expression estimates from all methods",
          metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## parameters specific to simulation experiments
# number of replications
n_rep <- opt$n_rep
i_rep <- opt$i_rep
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
if(cor_model %in% c('AR', 'MA', 'AR_10')){
	# correlation strength
	if(K == 2){
		rhos <- c(opt$rho1, opt$rho2)
	}else if(K == 4){
		rhos <- rep(0.8, 4)
	}else if(K == 10){
		rhos <- rep(0.8, 10)
	}
}
if(K == 2){
	# parameters for beta distribution, s.t. beta ~ beta(beta_1, beta_2)
	beta1 <- opt$beta1
	beta2 <- opt$beta2
	betas <- c(beta1, beta2)
}else if(K == 4){
	betas <- c(5, 2, 2, 1)
}else if(K == 10){
	betas <- c(4.5, 1.8, 1.8, 0.9, 0.3, 0.3, 0.1, 0.1, 0.1, 0.1)
}

## parameters for CSNet regressions
# methods <- list(var = opt$method_var,
# 	covar = opt$method_covar)

## parameters for sensitivity analysis
# whether to conduct sensitivity analysis
sensitivity <- opt$sensitivity == 'T'
if(sensitivity){
  print('Run sensitivity analysis')
  kappa <- opt$kappa
  b <- opt$b
}

## parameters specific to K=4 analysis
# whether to specify equal strength
equal_strength <- opt$equal_strength == 'T'

## whether to save co-expression estimates
save_est <- opt$save_est == 'T'

## fixed parameters

# set S=1 million s.t. simulated counts can be interpreted as FPKM
# (not considering the heterogeneity in gene lengths)
# set Tsize to be a smaller integer
# to decrease the level of bias in simulated correlations
NB_exper_pars <- list(Tsize=200, S=6e+07)

cv_mode <- 'independent'
cor_struct <- 'three_clusters'
gen_th_op <- 'SCAD'
n_val <- n #opt$n_val # number of samples for independent validation
exp_model <- 'NB'

print(log_var)
print(equal_strength)
# set file prefix
#if(K == 2 | K == 4){
	if(cor_model %in% c('AR', 'MA', 'AR_10')){
		if(K == 2){
			cor_model_prefix <- sprintf('%s_rho_%.1f_%.1f', cor_model, rhos[1], rhos[2])
		}else{
			cor_model_prefix <- cor_model
		}
	}else{
		cor_model_prefix <- cor_model
	}
	occ_prefix <- sprintf('log_var_%.1f_equal_strength_%s', log_var, equal_strength)
	prefix <- sprintf('%s/K_%i/%s/n_%i_p_%i', 
		cor_model_prefix, K,  
                occ_prefix,
		# ifelse((log_var != 8.0) | (!equal_strength), occ_prefix, 'standard'),
		n, p)
#}


if(sensitivity){
	prefix <- sprintf('%s/%s/kappa_%.2f_b_%.2f', 'sensitivity', prefix, kappa, b)
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

print(warnings())

# simulate expression data & run methods
# keep track of errors
coexp_methods <- c('d-Bulk', 'Bulk', 'd-CSNet-ols', 'CSNet-ols', 'd-CSNet', 'CSNet', 'bMIND-inf', 's-bMIND-inf', 'bMIND', 's-bMIND', 'ENIGMA', 's-ENIGMA', 'oracle')#, 'ENIGMA-trace', 's-ENIGMA-trace')
if(K == 10) coexp_methods <- coexp_methods[c(1:6, 13)]
metrics <- c('F_norm', 'op_norm', 'FPR', 'TPR')
error_mat <- matrix(NA, nrow = 4, ncol = K)
colnames(error_mat) <- paste0('cell_type_', 1:K)
rownames(error_mat) <- metrics

error_rec <- lapply(1:length(coexp_methods), function(i) error_mat)
names(error_rec) <- coexp_methods

est_list <- list()

# for(i_rep in 1:n_rep){
  # -
  # simulate data
  # -
  # simulate cell-type-specific expression
  train_list <- sim_exp(cor_model, n, seeds[i_rep], sim_setting, props = NULL, verbose=F)
  valid_list <- sim_exp(cor_model, n_val, seeds[i_rep]*2 + 10, sim_setting, 
    props = sim_setting$props[1:n_val, ],verbose=F)
 
  if(sensitivity){ 
    # update the proportions to be 'working' proportions (w/ bias & noise)
    wk_train_P <- add_noise_to_pi(train_list$data$P, b, kappa, 1)
    train_list$data$P <- wk_train_P
    valid_list$data$P <- wk_train_P[1:n_val, ]
    sim_setting$props <- wk_train_P
  }
  
  # evaluate ct-specific co-expressions estimated by different methods
  # -
  # evaluate CSNet
  # -
  
  methods_settings <- list(
    ols = list(var = 'ols', covar = 'ols'),
    wls = list(var = 'nnls', covar = 'wls')
  )
  settings_names <- c('-ols', '')
  
  for(i in 1:2){
    # obtain d-CSNet estimates
    train_cov_est <- mom_ls(train_list$data$P, train_list$data$X, methods_settings[[i]])
    train_cor_est <- lapply(train_cov_est, function(x) get_cor_from_cov(x)) 
    error_rec[[sprintf('d-CSNet%s', settings_names[i])]] <- sapply(1:K, function(k){
      eval_errors(train_cor_est[[k]], sim_setting$R[[k]], metrics, dense_estimate = T)
    })
    est_list[[sprintf('d-CSNet%s', settings_names[i])]] <- train_cor_est
    # tune thresholding parameter with the independent validation data
    tuning_result <- th_tuning(train_list$data, valid_list$data,
      coexp_method = 'CSNet',
      ctrl_list = list(ls_methods = methods_settings[[i]]))
  
    # obtain CSNet estiamtes
    train_cor_th_est <- lapply(1:K, function(k){
      generalized_th(train_cor_est[[k]], tuning_result$th[k], gen_th_op, F)
    })
    error_rec[[sprintf('CSNet%s', settings_names[i])]] <- sapply(1:K, function(k){
      eval_errors(train_cor_th_est[[k]], sim_setting$R[[k]], metrics)
    })
    est_list[[sprintf('CSNet%s', settings_names[i])]] <- train_cor_th_est
  }
  
  # -
  # evaluate sparse bulk estimate
  # -
  bulk_train_cor_est <- cor(train_list$data$X)
  bulk_tuning_result <- th_tuning(train_list$data, valid_list$data, 
    coexp_method = 'Bulk')
  bulk_train_cor_th_est <- lapply(1:K, function(k){
    generalized_th(bulk_train_cor_est, bulk_tuning_result$th, gen_th_op, F)
  }) # same for two cell types
  error_rec[['Bulk']] <- sapply(1:K, function(k){
    eval_errors(bulk_train_cor_th_est[[k]], sim_setting$R[[k]], metrics)
  })
  est_list[['Bulk']] <- bulk_train_cor_th_est
  est_list[['d-Bulk']] <- lapply(1:K, function(k) bulk_train_cor_est)
  
  # -
  # evaluate oracle estimate
  # -
  oracle_train_cor_est <- lapply(1:K, function(k) cor(train_list$ct_specific_data[[k]]))
  oracle_tuning_result <- lapply(1:K, function(k){
    th_tuning(list(X=train_list$ct_specific_data[[k]]), list(X=valid_list$ct_specific_data[[k]]), coexp_method = 'Bulk')})
  oracle_train_cor_th_est <- lapply(1:K, function(k){
    generalized_th(oracle_train_cor_est[[k]], oracle_tuning_result[[k]]$th, gen_th_op, F)
  })
  error_rec[['d-oracle']] <- sapply(1:K, function(k){
    eval_errors(oracle_train_cor_est[[k]], sim_setting$R[[k]], metrics, dense_estimate = T)
  })
  error_rec[['oracle']] <- sapply(1:K, function(k){
    eval_errors(oracle_train_cor_th_est[[k]], sim_setting$R[[k]], metrics)
  })
  est_list[['d-oracle']] <- oracle_train_cor_est
  est_list[['oracle']] <- oracle_train_cor_th_est

  print(warnings())

  # -
  # evaluate bMIND
  # -
  if(K != 10){
  ## with informative prior
  prior_info_vec <- c('informative', 'non-informative')
  suffix_vec <- c('-inf', '')
  
  for(i in 1:2){
    # bMIND
    bMIND_train_est <- coexp_by_bMIND(train_list, prior_info = prior_info_vec[i])
    error_rec[[sprintf('bMIND%s', suffix_vec[i])]] <- sapply(1:K, function(k){
      eval_errors(bMIND_train_est[[k]], sim_setting$R[[k]], metrics, dense_estimate = T)
    })
    est_list[[sprintf('bMIND%s', suffix_vec[i])]] <- bMIND_train_est
    
    if(anyNA(error_rec[[sprintf('bMIND%s', suffix_vec[i])]])){
      error_rec[[sprintf('s-bMIND%s', suffix_vec[i])]] <- matrix(NA, nrow = length(metrics), ncol = K)
      est_list[[sprintf('s-bMIND%s', suffix_vec[i])]] <- lapply(1:K, function(k) diag(1, p))
    }else{
      # sparse bMIND
      bMIND_tuning_result <- th_tuning(train_list, valid_list,
      coexp_method = 'bMIND',
      ctrl_list = list(prior = prior_info_vec[i]))
      sparse_bMIND_train_th_est <- lapply(1:K, function(k){
        generalized_th(bMIND_train_est[[k]], bMIND_tuning_result$th[k], gen_th_op, F)
      })
      error_rec[[sprintf('s-bMIND%s', suffix_vec[i])]] <- sapply(1:K, function(k){
        eval_errors(sparse_bMIND_train_th_est[[k]], sim_setting$R[[k]], metrics)
      })
      est_list[[sprintf('s-bMIND%s', suffix_vec[i])]] <- sparse_bMIND_train_th_est
    }
  }
  
  # -
  # evaluate ENIGMA
  # -
  ## with L2 max norm
  # ENIGMA
  ENIGMA_train_est <- coexp_by_ENIGMA(train_list, norm = 'L2', seed = seeds[i_rep])
  error_rec[['ENIGMA']] <- sapply(1:K, function(k){
    eval_errors(ENIGMA_train_est[[k]], sim_setting$R[[k]], metrics, dense_estimate = T)
  })
  est_list[['ENIGMA']] <- ENIGMA_train_est  

  # sparse ENIGMA
  ENIGMA_tuning_result <- th_tuning(train_list, valid_list,
    coexp_method = 'ENIGMA',
    ctrl_list = list(norm = 'L2', seed = seeds[i_rep]))
  sparse_ENIGMA_train_th_est <- lapply(1:K, function(k){
    generalized_th(ENIGMA_train_est[[k]], ENIGMA_tuning_result$th[k], gen_th_op, F)
  })
  error_rec[['s-ENIGMA']] <- sapply(1:K, function(k){
    eval_errors(sparse_ENIGMA_train_th_est[[k]], sim_setting$R[[k]], metrics, dense_estimate = F)
  })
  est_list[['s-ENIGMA']] <- sparse_ENIGMA_train_th_est
  
  ## with trace norm
  # ENIGMA
  #ENIGMA_trace_train_est <- coexp_by_ENIGMA(train_list, norm = 'trace', seed = seeds[i_rep])
  #error_rec[['ENIGMA-trace']] <- sapply(1:K, function(k){
  #  eval_errors(ENIGMA_trace_train_est[[k]], sim_setting$R[[k]], metrics, dense_estimate = T)
  #})
  #est_list[['ENIGMA-trace']] <- ENIGMA_trace_train_est

  # sparse ENIGMA
  #ENIGMA_trace_tuning_result <- th_tuning(train_list, valid_list,
  #  coexp_method = 'ENIGMA',
  #  ctrl_list = list(norm = 'trace', seed = seeds[i_rep]))
  #sparse_ENIGMA_trace_train_th_est <- lapply(1:K, function(k){
  #  generalized_th(ENIGMA_trace_train_est[[k]], ENIGMA_trace_tuning_result$th[k], gen_th_op, F)
  #})
  #error_rec[['s-ENIGMA-trace']] <- sapply(1:K, function(k){
  #  eval_errors(sparse_ENIGMA_trace_train_th_est[[k]], sim_setting$R[[k]], metrics, dense_estimate = T)
  #})
  #est_list[['s-ENIGMA-trace']] <- sparse_ENIGMA_trace_train_th_est

# }
  }

saveRDS(error_rec, sprintf('%s/n_rep_%i_i_rep_%i_error.rds', result_prefix, n_rep, i_rep))
if(save_est){
  saveRDS(est_list, sprintf('%s/n_rep_%i_i_rep_%i_est.rds', result_prefix, n_rep, i_rep))
  if(i_rep==1) saveRDS(sim_setting, sprintf('%s/sim_setting.rds', result_prefix))
}

if(n_rep == 1){
  g_list <- list()
  t <- 1
  for(coexp_method in coexp_methods){
    for(k in 1:K){
	if(K != 10){
		g <- plot_heatmap(est_list[[coexp_method]][[k]], sim_setting$sub_cl, 
             title = sprintf('%s, ct %i', coexp_method, k))
	}else{
		g <- plot_heatmap(est_list[[coexp_method]][[k]], 
				  lapply(1:4, function(j) sim_setting$sub_cl[[j]]),
				  title = sprintf('%s, ct %i', coexp_method, k))
	}
	g_list[[t]] <- g[[4]]
        ggsave(sprintf('%s/%s_ct_%i.pdf', fig_prefix, coexp_method, k), g[[4]])
        t <- t + 1
    }
  }
  saveRDS(g_list, sprintf('%s/all_methods.rds', fig_prefix))
}
#ggsave(sprintf('%s/all_methods.pdf', fig_prefix),
#       grid.arrange(grobs = g_list, nrow = length(coexp_methods), ncol = K),
#       width = 7 * K, height = 7 * length(coexp_methods))
