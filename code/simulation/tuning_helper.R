# 
# tuning_helper.R
# 
# tune the thresholding parameters in simulation
# with independently generated validation data

th_tuning <- function(train_data, valid_data,
	coexp_method,
	ctrl_list = NULL, th_func = 'SCAD', n_th = 100, to_plot = F){
	source('../helper_function.R')
	# evaluate estimates on training and validation data
	if(coexp_method == 'CSNet'){
		train_est <- mom_ls(train_data$P, train_data$X, ctrl_list$ls_methods)
		valid_est <- mom_ls(valid_data$P, valid_data$X, ctrl_list$ls_methods)
		K <- ncol(train_data$P)
		# threshold correlation matrices instead of covariance matrices
		train_est <- lapply(train_est, function(x) get_cor_from_cov(x))
		valid_est <- lapply(valid_est, function(x) get_cor_from_cov(x))
	}else if(coexp_method == 'Bulk'){
		train_est <- cor(train_data$X) %>% list
		valid_est <- cor(valid_data$X) %>% list
		K <- 1
	}else if(coexp_method == 'bMIND'){
		train_est <- coexp_by_bMIND(train_data, prior_info = ctrl_list$prior)
		valid_est <- coexp_by_bMIND(valid_data, prior_info = ctrl_list$prior)
		K <- ncol(train_data$data$P)
	}else if(coexp_method == 'ENIGMA'){
		train_est <- coexp_by_ENIGMA(train_data, norm = ctrl_list$norm, seed = ctrl_list$seed)
		valid_est <- coexp_by_ENIGMA(valid_data, norm = ctrl_list$norm, seed = ctrl_list$seed*2 + 19)
		K <- ncol(train_data$data$P)
	}
	
	# tune thresholding parameters 
	# by minimizing F norm
	F_norm_rec <- array(0, c(K, n_th, 2))
	th_grids <- sapply(1:K, function(k) seq(0, 1, length.out = n_th))
	for(k in 1:K){
		for(j in 1:n_th){
			train_est_k <- generalized_th(train_est[[k]], th_grids[j, k], 
				th_func, th_diag = F)
			F_norm_rec[k, j, 1] <- norm(valid_est[[k]] - train_est_k, 'F')
		}
	}
	# Found a bug
	# F_norm_rec should have [,,2] from thresholding valid est
	F_norm_means <- apply(F_norm_rec, c(1,2), function(x) mean(x))
	F_norm_sds <- apply(F_norm_rec, c(1,2), function(x) sd(x))

	# select the threshold which minimizes average CV error
	selected_th <- numeric(K)
	for(k in 1:K){
    	selected_th[k] <- th_grids[which.min(F_norm_means[k,]) , k]
    }

    if(coexp_method != 'Bulk'){ 
    	ct_names <- colnames(train_data$P)
    }else{
    	ct_names <- 'bulk'
    }

    if(to_plot){
        tmp_df <- data.frame(th = as.vector(th_grids),
                         ct = rep(ct_names, each = nrow(th_grids)),
                         F_norm = as.vector(t(F_norm_means)),
                         sd = as.vector(t(F_norm_sds)),
                         selected_th = rep(selected_th, each = nrow(th_grids)))
        g <- ggplot(tmp_df, aes(x = th, y = F_norm)) +
        geom_point() +
        geom_line(alpha = 0.8) +
        geom_errorbar(aes(ymin = F_norm-sd, ymax = F_norm+sd), alpha = 0.5) +
        geom_vline(data=tmp_df, aes(xintercept=selected_th)) +
        facet_wrap(~ct) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

        print(g)
    }
    return(list(F_norm_means = F_norm_means,
    	F_norm_sds = F_norm_sds,
    	th_grids = th_grids,
    	th = selected_th))
}
