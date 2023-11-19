##
# sample level perturbation
## 
sample_level_perturb <- function(kappa, P, seed=1){
    require(DirichletReg)	
    set.seed(seed)
    n_train <- nrow(P)
    props_0 <- matrix(0, nrow = n_train, ncol = 7)
    for(i in 1:n_train){
        pi_i <-  P[i,]
        if(any(pi_i == 0)){
            # impute zero proportion with a small constant
            # as alpha can not have zero component for Dirichlet distribution
            zero_ind <- which(pi_i == 0)
            pi_i[zero_ind] <- 2e-3
        }
        props_0[i,] <- rdirichlet(1, pi_i * kappa)
    }
    colnames(props_0) <- colnames(P)
    rownames(props_0) <- rownames(P)
    return(props_0)
}

##
# systematic experiments
##
real_data_sens <- function(geneset_list, kappa_range = c(10), n_iter=2){
    # print(n_iter)
    n_kappa <- length(kappa_range)
    F_rec <- TPR_rec <- FPR_rec <- array(0, dim = c(7, n_kappa, n_iter))

    # print(dim(F_rec))

    est_list_list <- list()
    for(i in 1:n_kappa){
        kappa <- kappa_range[i]
        for(i_iter in 1:n_iter){
            # print(i_iter)
            seed = 87+i_iter**2

            # -
            # new: sample level pertubation
            # -
            props <- sample_level_perturb(kappa, geneset_list$full_data_list[[1]]$P, seed)
            # run our method
            methods_settings <- list(var='nnls', covar='wls')
            data_list <- geneset_list$full_data_list[[1]]
            dCSNet_cov_est <- mom_ls(props, data_list$X, methods_settings, weight_th = 0.01)
            dCSNet_est <- lapply(dCSNet_cov_est, function(x) get_cor_from_cov(x))
            cv_res <- cross_validation(data_list,
                           coexp_method = 'CSNet',
                           list(methods_settings = methods_settings, th = 0.01),
                           ct_names = colnames(data_list$P),
                           to_plot = F,
                          seed =  1)
            CSNet_est <- lapply(1:length(dCSNet_est), function(i){
                tmp <- dCSNet_est[[i]]
                na_inds <- is.na(tmp)
                tmp[na_inds] <- 0

                th_tmp <- generalized_th(tmp, cv_res$th[i], 'SCAD', F)
                th_tmp[na_inds] <- NA
                diag(th_tmp) <- 1
                return(th_tmp)
            })


            if(i_iter == 1){
                est_list_list[[as.character(kappa)]] <- CSNet_est
            }

            # evaluate difference with original estimates
            for(k in 1:7){
                F_rec[k, i, i_iter] <- F_norm(CSNet_est[[k]],
                                              geneset_list$CSNet_est[[k]])
                TPR_rec[k, i, i_iter] <- eval_TPR(CSNet_est[[k]],
                                                  geneset_list$CSNet_est[[k]],
                                                 off_diagonal = F)
                FPR_rec[k, i, i_iter] <- eval_FPR(CSNet_est[[k]],
                                                  geneset_list$CSNet_est[[k]],
                                                 off_diagonal = F)
            }

        }
        rownames(F_rec) <- rownames(TPR_rec) <- rownames(FPR_rec) <- colnames(geneset_list$full_data_list[[1]]$P)
    }
    return(list(F_norm = F_rec, TPR = TPR_rec, FPR = FPR_rec, est_list_list))
}

# for visualization and extract legend from a plot
library(grid)
library(gridExtra)
## Function to extract legend
g_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}

##
# plot metrics against kappa
##
vis_metric <- function(rec, metric, legend, n_kappa = length(kappa_range)){
    require(reshape2)
    rec_df <- melt(rec[1:4,1,])
    n_iter <- dim(rec)[3]
    rec_df$kappa <- rep((kappa_range)[1], 4*n_iter)
    for(j in 2:n_kappa){
        rec_tmp <- melt(rec[1:4,j,])
        rec_tmp$kappa <- rep((kappa_range)[j], 4*n_iter)
        rec_df <- rbind(rec_df, rec_tmp)
    }
    rec_df$kappa <- factor(rec_df$kappa)
    # print(head(rec_df))

    legend_pos <- ifelse(legend, 'right', 'none')

    g<-ggplot(rec_df, aes(x = kappa, y = value)) +
    geom_boxplot(aes(color = Var1, group = interaction(Var1, kappa))) +
    #geom_line(aes(group = Var2, color = Var1)) +
    labs(y = gsub('_', ' ', metric), title = gsub('_', ' ', metric),
         color = 'Cell type') +
    theme(text = element_text(size=13),
     plot.title = element_text(size=15),
         legend.position = legend_pos)
    if(metric %in% c('TPR', 'FPR')){
        g + ylim(c(0,1))
    }
    return(g)
}

##
# save senstivity result figures
##
make_fig <- function(sens_res, name, title){
    sens_g <- list()
    names(sens_res) <- c('F-norm', 'TPR', 'FPR', 'est')
    for(metric in c('F-norm', 'TPR', 'FPR')){
        sens_g[[metric]] <- vis_metric(sens_res[[metric]], metric, legend=F)
    }
    tmp <-  vis_metric(sens_res[[metric]], metric, legend=T)
    sens_g[['legend']] <- g_legend(tmp)
    if(is.null(title)){
    combined_g <- grid.arrange(grobs = sens_g, ncol=4, nrow=1,
                             widths = c(3,3,3,1))
    }else{
    combined_g <- grid.arrange(grobs = sens_g, ncol=4, nrow=1,
                             widths = c(3,3,3,1),
                             top = textGrob(sprintf("%s",
                                                   #gsub('_', ' ', tolower(name))
                                                   title),
                                            gp=gpar(fontsize=22,font=3)))
}

    ggsave(sprintf('%s/sensitivity/%s_robustness_metrics.pdf', 'figures', name),
       combined_g,
       width = 10, height = 3)
}

