
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

cross_validation <- function(data_list, methods_settings = NULL, cv_fold = 2,
                             n_th = 100, th_func = 'SCAD',
                             seed = 43978, to_plot = F, ct_names = NULL){
    ct_est <- ifelse(length(data_list) == 1, F, T)
    K <- ifelse(ct_est, ncol(data_list$P), 1)
    n <- nrow(data_list$X)
    th_grids <- sapply(1:K, function(k) seq(0, 1, length.out = n_th))
    set.seed(seed)
    perm_inds <- sample.int(n, n)
    assign_inds <- cut(seq(1,n), breaks=cv_fold, labels=FALSE)
    F_norm_rec <- array(0, dim = c(K, n_th, cv_fold))
    for(i in 1:cv_fold){
        val_inds <- perm_inds[assign_inds == i]
        train_inds <- perm_inds[assign_inds != i]
        if(ct_est){
            train_est <- mom_ls(data_list$P[train_inds,], data_list$X[train_inds,], methods_settings)
            val_est <- mom_ls(data_list$P[val_inds, ], data_list$X[val_inds, ], methods_settings)
            val_est <- lapply(val_est, function(est) get_cor_from_cov(est))
            train_est <- lapply(train_est, function(est) get_cor_from_cov(est))
        }else{
            train_est <- cor(data_list$X[train_inds,] %>% t) %>% list()
            val_est <- cor(data_list$X[val_inds, ] %>% t) %>% list()
        }
        # tune thresholding parameters
        for(k in 1:K){
        for(j in 1:n_th){
          train_est_k <- generalized_th(train_est[[k]], th_grids[j, k], th_func, th_diag = F)
          val_est_k <- val_est[[k]]
          F_norm_rec[k, j, i] <- norm(val_est_k - train_est_k,'F')
        }
      }
    }
    F_norm_means <- apply(F_norm_rec, c(1,2), function(x) mean(x))
    F_norm_sds <- apply(F_norm_rec, c(1,2), function(x) sd(x))
    # select the threshold which minimizes average CV error
    selected_th <- numeric(K)
    for(k in 1:K){
      selected_th[k] <- th_grids[which.min(F_norm_means[k,]) , k]
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
