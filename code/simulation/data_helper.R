# 
# simulate_exp.R
# 
# simualate cell-type-specific expressions and assemble them into bulk expressions

sim_exp <- function(cor_model, n, seed, sim_setting, props = NULL, verbose = F){
  source('NB_helper.R')
  # extract simulation parameters
  p <- length(sim_setting[['mu']][[1]])
  NB_par_list <- sim_setting[['NB_par']]
  mu_star_list <- sim_setting[['mu']]
  sigma_sq_star_list <- sim_setting[['sigma_sq']]
  R_star_list <- sim_setting[['R']]
  K <- length(sim_setting[['mu']])
  cor_model <- sim_setting[['cor_model']]
  # dimension and ids for cell-type-specific gene sets
  cor_p <- sim_setting[['cor_p']]
  sub_cl <- sim_setting[['sub_cl']]

  gene_exp_list <- list()
  for(k in 1:K){
    # make sure the random seed for two cell types are different
    ct_specific_seed <- seed + (k^2+43)
    gene_exp_mat <- matrix(0, n, p)
    
    # generate correlated expressions
    if(cor_model != 'dcSBM'){
    	tmp_cor <- simu_NB_obs(n,
                           NB_par_list[[k]],
                           1:cor_p,
                           model=cor_model,
                           seed=ct_specific_seed)
    }
    gene_exp_mat[, sub_cl[[k]]] <- tmp_cor

    # generate independent NB observations
    tmp_indpt <- matrix(NA, n, p-cor_p)
    # extract mean and variance statistics for independent genes
    indpt_ind <- which(!(1:p %in% sub_cl[[k]]))
    gene_mu_sub <- mu_star_list[[k]][indpt_ind]
    gene_var_sub <- sigma_sq_star_list[[k]][indpt_ind]
    set.seed(ct_specific_seed)
    for(j in 1:(p-cor_p)){
      gene_mu <- gene_mu_sub[j]
      gene_var <- gene_var_sub[j]
      gene_p <- gene_mu / gene_var
      tmp_indpt[,j] <- stats::rnbinom(n, 
        size = gene_mu * gene_p / (1-gene_p),
        mu = gene_mu)
    }
    gene_exp_mat[, indpt_ind] <- tmp_indpt
    gene_exp_list[[k]] <- gene_exp_mat
    }
    # return(gene_exp_list)
    # obtain bulk expressions
    bulk_exp <- matrix(0, nrow = n, ncol = p)
    for(k in 1:K){
      bulk_exp <- bulk_exp + gene_exp_list[[k]] * sim_setting$props[,k]
    }
    if(is.null(props)) props <- sim_setting$props
    return(list(data = list(X = bulk_exp, P = props), ct_specific_data = gene_exp_list))
}