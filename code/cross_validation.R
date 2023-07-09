
# visualize cv results
vis_cv <- function(cv_res){
  cv_g_list <- list()
  for(k in 1:2){
    cv_df <-data.frame(mean = cv_res$F_norm_means[k,], sd = cv_res$F_norm_sds[k,],
                       th = cv_res$th_grids[,k])
    cv_g_list[[k]] <- ggplot(cv_df, aes(x = th, y = mean)) +
      geom_point() +
      geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd),  alpha = 0.2, width = 0.01) +
      labs(title = k, y = 'F norm')
  }

  cv_g <- grid.arrange(grobs = cv_g_list, nrow = 1, top = 'CV error curve')
  return(cv_g)
}

cross_validation <- function(data_list, 
			     coexp_method = 'CSNet',
			     method_pars = NULL,
			     cv_fold = 2,
                             n_th = 100, 
			     th_func = 'SCAD',
                             seed = 1, 
			     to_plot = F, 
			     ct_names = NULL,
			     n_splits = 1){
    ct_est <- ifelse(length(data_list) == 1, F, T)
    K <- ifelse(ct_est, ncol(data_list$P), 1)
    n <- nrow(data_list$X)
    th_grids <- sapply(1:K, function(k) seq(0, 1, length.out = n_th))
    set.seed(seed)
    
    print(sprintf('n splits = %i', n_splits))
    F_norm_rec_large <-  array(0, dim = c(K, n_th, n_splits))

    for(i_split in 1:n_splits){

    perm_inds <- sample.int(n, n)
    assign_inds <- cut(seq(1,n), breaks=cv_fold, labels=FALSE)
    F_norm_rec <- array(0, dim = c(K, n_th, cv_fold))
    for(i in 1:cv_fold){
        val_inds <- perm_inds[assign_inds == i]
        train_inds <- perm_inds[assign_inds != i]
        if(ct_est){
		if(coexp_method == 'CSNet'){
        		train_est <- mom_ls(data_list$P[train_inds,], data_list$X[train_inds,], 
					    method_pars$methods_settings,
					    method_pars$th)
        		val_est <- mom_ls(data_list$P[val_inds, ], data_list$X[val_inds, ], 
					  method_pars$methods_settings, 
					  method_pars$th)
        		val_est <- lapply(val_est, function(est) get_cor_from_cov(est))
        		train_est <- lapply(train_est, function(est) get_cor_from_cov(est))
		}else if(coexp_method == 'bMIND'){
			# source()
			train_list <- list(X = data_list$X[train_inds,], P = data_list$P[train_inds, ])
			train_est <- real_data_run_bMIND(data_list = train_list, 
							 method_pars$ind_bulk, method_pars$ind_profile)
			val_list <- list(X = data_list$X[val_inds, ], P = data_list$P[val_inds, ])
			val_est <- real_data_run_bMIND(data_list = val_list,
							 method_pars$ind_bulk, method_pars$ind_profile)
		}
	}else{
            train_est <- cor(data_list$X[train_inds,]) %>% list()
            val_est <- cor(data_list$X[val_inds, ]) %>% list()
        }
        # tune thresholding parameters
        for(k in 1:K){
        for(j in 1:n_th){
	    # exclude NA entries when evaluating F norm
	    train_non_NA <- train_est[[k]]
	    train_non_NA[is.na(train_non_NA)] <- 0
	    diag(train_non_NA) <- 1
            train_est_k <- generalized_th(train_non_NA, th_grids[j, k], th_func, th_diag = F)
            val_est_k <- val_est[[k]]
	    val_est_k[is.na(val_est_k)] <- 0
	    diag(val_est_k) <- 1
            F_norm_rec[k, j, i] <- norm(val_est_k - train_est_k,'F')
        }
      }
    }
    F_norm_rec_large[, , i_split] <- F_norm_rec[,,1] + F_norm_rec[,,2]
    }

    F_norm_means <- apply(F_norm_rec_large, c(1,2), function(x) mean(x))
    F_norm_sds <- apply(F_norm_rec_large, c(1,2), function(x) sd(x))
    F_norm_ses <- F_norm_sds / sqrt(n_splits)
    # select the threshold which minimizes average CV error
    selected_th <- selected_th_1_se <- numeric(K)
    for(k in 1:K){
	min_ind <- which.min(F_norm_means[k,])
        selected_th[k] <- th_grids[min_ind , k]
      
        if(n_splits > 1){
	    one_se_ind <- min(which(F_norm_means[k,] <= F_norm_means[k, min_ind] + F_norm_ses[k, min_ind]))
	    selected_th_1_se[k] <- th_grids[one_se_ind, k]
        }
    }
    if(K > 1 & is.null(ct_names)){
        ct_names <- paste0('cell type ', 1:K)
    }
    if(to_plot){
    # plot CV curves
    tmp_df <- data.frame(th = as.vector(th_grids),
                         ct = rep(ct_names, each = nrow(th_grids)),
                         F_norm = as.vector(t(F_norm_means)),
                         sd = as.vector(t(F_norm_sds)),
                         selected_th = rep(selected_th, each = nrow(th_grids)),
			 selected_1se_th = rep(selected_th_1_se, each = nrow(th_grids)))
    g <- ggplot(tmp_df, aes(x = th, y = F_norm)) +
      geom_point() +
      geom_line(alpha = 0.8) +
      geom_errorbar(aes(ymin = F_norm-sd, ymax = F_norm+sd), alpha = 0.5) +
      geom_vline(data=tmp_df, aes(xintercept=selected_th)) +
      # geom_vline(data=tmp_df, aes(xintercept=selected_1se_th), color = 'blue') +
      facet_wrap(~ct) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    if(n_splits > 1) g <- g + geom_vline(data=tmp_df, aes(xintercept=selected_1se_th), color = 'blue')

    print(g)
    }

    return(list(F_norm_means = F_norm_means,
              F_norm_sds = F_norm_sds,
              th_grids = th_grids,
              th = selected_th,
	      th_1se = selected_th_1_se,
	      train_est = train_est,
	      val_est = val_est))
}
