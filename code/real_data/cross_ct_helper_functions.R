# compute pseudo-bulk data based on single cell data
load_data <- function(ct){
    # counts
    counts <- readMM(sprintf('%s/%s_counts.txt', data_dir, ct))
    features <- read.table(sprintf('%s/%s_features.txt', data_dir, ct))[[1]]
    subjects <- read.table(sprintf('%s/%s_subjects.txt', data_dir, ct))[[1]]
    rownames(counts) <- features
    colnames(counts) <- subjects
    # sparisty
    sparsity <- Matrix::rowSums(counts == 0) / ncol(counts)
    # seq depth
    seq_depth <- Matrix::colSums(counts)
    # remove samples with no reads
    counts <- counts[, - which(seq_depth == 0)] # same for all four cell types
    seq_depth <- seq_depth[- which(seq_depth == 0)]
    # log normalized data
    dense_counts <- as.matrix(counts)
    scaled_data <- t(t(dense_counts) / seq_depth) * median(seq_depth)
    log_nor_data <- log(scaled_data + 1)
    # mean expression levels
    mean_exp <- rowMeans(log_nor_data)
    mean_exp_order <- length(mean_exp) - rank(mean_exp)
    return(list(counts = counts, log_nor_data = log_nor_data, scaled_data = scaled_data,
               sample_summary = list(seq_depth = seq_depth),
               gene_summary = list(sparsity = sparsity, mean_exp = mean_exp, mean_exp_order = mean_exp_order)))
}

load_genes <- function(geneset){
    # load relevant data
    if(geneset != 'AD_risk_genes'){
        genes <- readRDS(sprintf('%s/%s_genes.rds', output_prefix, geneset))
    }else{
        genes <- readRDS(sprintf('%s/AD_genes_anlyzed.rds', '~/project/CSNet/real_data/covest-real-data/ROSMAP/output/AD_genes/multiplex_new'))
    }
    return(genes)
}

# compare the magnitude of within and between
compare_cross_and_within <- function(pseudo_bulk_sc, weighting = F, filtering = F, method = 'pearson',
                                    gaussian_matching = F){
    K <- length(pseudo_bulk_sc)
    p <- nrow(pseudo_bulk_sc[[1]])
    within_sum <- between_sum <- matrix(0, p, p)
    cor_p_val <- cov_all_cts <- array(NA, dim = c(p, p, K, K))
    ct_props_matched <- ct_props[match(names(pseudo_bulk_sc), names(ct_props))]
    print(ct_props_matched)
    if(gaussian_matching){
        for(k1 in 1:K){
            for(i in 1:p){
                tmp <- pseudo_bulk_sc[[k1]][i, ]
                tmp_qt <- rank(tmp) / length(tmp)
                tmp_qm <- qnorm(tmp_qt, mean(tmp), sd(tmp))
                tmp_qm[is.infinite(tmp_qm)] <- mean(tmp) + 3* sd(tmp)
                pseudo_bulk_sc[[k1]][i, ] <- tmp_qm

            }
        }
    }
    # gene 1
    for(i in 1:p){
        # gene 2
        for(j in 1:p){
            # cell type 1
            for(k1 in 1:K){
                # within-cell-type
                exp1 <- pseudo_bulk_sc[[k1]][i, ]
                exp2 <- pseudo_bulk_sc[[k1]][j, ]
                if(sd(exp1) > 0 & sd(exp2) > 0){
                    if(weighting){
                        cov_12 <- cov(exp1, exp2) * ct_props_matched[k1]^2
                    }else{
                        cov_12 <- cov(exp1, exp2)
                    }
                    pval <- cor.test(exp1, exp2, method = method)$p.value
                    if(filtering){
                        cov_12 <- ifelse(pval < 0.05, cov_12, 0)
                    }
                    cov_all_cts[i, j, k1, k1] <- cov_12
                    within_sum[i, j] <- within_sum[i, j] + cov_12
                    cor_p_val[i, j, k1, k1] <- pval
                }

                for(k2 in 1:K){
                    # between-cell-type covariance
                    if(k1 != k2){
                        exp1 <- pseudo_bulk_sc[[k1]][i, ]
                        exp2 <- pseudo_bulk_sc[[k2]][j, ]
                        if(sd(exp1) > 0 & sd(exp2) > 0){
                            if(weighting){
                                cov_12 <- cov(exp1, exp2) * ct_props_matched[k1]  * ct_props_matched[k2]
                            }else{
                                cov_12 <- cov(exp1, exp2)
                            }
                            pval <- cor.test(exp1, exp2, method = method)$p.value
                            if(filtering){
                                cov_12 <- ifelse(pval < 0.05, cov_12, 0)
                            }
                            between_sum[i, j] <- between_sum[i, j] + cov_12
                            cov_all_cts[i, j, k1, k2] <- cov_12
                            cor_p_val[i, j, k1, k2] <- pval
                        }
                    }
                }
            }
        }
    }
    return(list(within = within_sum, between = between_sum, cov_rec = cov_all_cts,
               cor_p_val = cor_p_val))
}

extract_within_and_between <- function(within_cov){
    within_ct_ind <- array(F, dim = dim(within_cov))
between_ct_ind <- array(T, dim = dim(within_cov))
for(i in 1:dim(within_cov)[1]){
    for(j in 1:dim(within_cov)[2]){
        for(k1 in 1:dim(within_cov)[3]){
            for(k2 in 1:dim(within_cov)[4]){
                if(k1 == k2){
                    within_ct_ind[i, j, k1, k2] <- T
                    between_ct_ind[i, j, k1, k2] <- F
                }
            }
        }
    }
}
    return(list(within = within_ct_ind, between = between_ct_ind))
}

# for making plots

eval_within_and_between <- function(res, inds, method = 'BH', p_cutoff = 0.05, thresholding = T, geneset = '',
                                   min_threshold = 1e-6, title = '', to_save=F){
    if(method == 'BH'){
        within_p_adjusted <- p.adjust(res$cor_p_val, method = 'BH')[inds$within]
        between_p_adjusted <- p.adjust(res$cor_p_val, method = 'BH')[inds$between]
    }else{
        within_p_adjusted <- res$cor_p_val[inds$within]
        between_p_adjusted <- res$cor_p_val[inds$between]
    }
    ngenes <- dim(inds$within)[1]
    p_cutoff <- ifelse(method == 'BH', p_cutoff, p_cutoff / (ngenes^2))
    if(!thresholding) p_cutoff = 1
    within <- res$cov_rec
    within[!inds$within] <- 0
    within[inds$within][within_p_adjusted > p_cutoff] <- 0
    within_sum <- apply(within, c(1,2), function(x) sum(x,na.rm=T))

    between <- res$cov_rec
    between[!inds$between] <- 0
    between[inds$between][between_p_adjusted > p_cutoff] <- 0
    between_sum <- apply(between, c(1,2), function(x) sum(x,na.rm=T))


    g2 <- #ggplot(cov_sum_df[(abs(cov_sum_df$within) > min_threshold & abs(cov_sum_df$between) > min_threshold), ]) +
    ggplot(cov_sum_df[(cov_sum_df$within !=0 & cov_sum_df$between != 0), ]) +
      geom_histogram(aes(x = log10(within/between)), alpha = 0.2, bins = 30) +
                         geom_vline(xintercept = 0, color = 'grey', linetype = 'dashed') +
                         geom_vline(aes(xintercept = median(log10(within / between) )),
                                    color = 'blue', linetype = 'dashed') +
                         labs(x = 'log10(|within|/|between|)') +
                         theme_classic(base_size = 20)
    print(g2)
    print(sprintf('median:%.2f', median(log10(cov_sum_df$within / cov_sum_df$between))))
          ggsave(sprintf('figures/%s_cross_vs_within%s.pdf', geneset, ifelse(thresholding, '_thresholding', '')),
                 g2 + labs(title = gsub('_', ' ', geneset)) + theme_classic(base_size = 20) +
                         theme(plot.title = element_text(hjust = 0.5)) , width = 6, height = 4)


    sprintf('method: %s, prop of within equal to between: %.2f; & greater than between: %.2f (%.2f)',
           method, mean((within_sum %>% abs) == (between_sum %>% abs), na.rm=T),
            mean((within_sum %>% abs) > (between_sum %>% abs), na.rm=T),
           sum((within_sum %>% abs) > (between_sum %>% abs), na.rm=T) / (sum((within_sum %>% abs) != (between_sum %>% abs), na.rm=T))) %>% print
    # return(cov_sum_df)
}


plot_within_versus_between <- function(gene_set, filtering, gene_set_name){
    if(filtering){
        add_inds <- top_inds[['10000']]
    }else{
        add_inds <- rep(T, length(features))
    }
    scaled_data_list <- lapply(major_cts, function(ct){
        data_list[[ct]]$scaled_data[features %in% gene_set & add_inds, ]
    })
    names(scaled_data_list) <- major_cts

    print(sprintf('Total genes: %i, after intersection with genes highly expressed in all four cell types: %i',
                  length(gene_set), nrow(scaled_data_list[[1]])))

    res <- compare_cross_and_within(scaled_data_list, weighting = T, filtering = F, method = 'pearson',
                                  gaussian_matching = F)
    inds <- extract_within_and_between(res$cov_rec)

    eval_within_and_between(res, inds, thresholding=F, title = sprintf('%s genes, no thresholding', gene_set_name),
                            geneset = gene_set_name, to_save = ifelse(filtering, T, F))
}

# Implementation of cross-CSNet
mom_ls_with_cross <- function (P, X, methods = list(var = "nnls", covar = "wls")) 
{
    p <- ncol(X)
    K <- ncol(P)
    n <- nrow(X)
    # P_2 <- P^2
    P_2 <- matrix(0, nrow = n, ncol = K*(K+1)/2)
    t <- 1
    for(i in 1:K){
        for(j in i:K){
            P_2[,t] <- P[,i] * P[,j]
            if(i == j) print(t)
            t <- t + 1
        }
    }
    mu <- apply(X, 2, function(x) nnls(P, x)$x)
    M <- P %*% mu
    Sigma_array <- array(NA, c(p, p, K*(K+1)/2))
    X_centered <- X - M
    Y <- X_centered^2
    if (methods$var == "nnls") {
        #sigma_var <- 
        #for(i in c(1,5,8,10))
        sigma_var <- apply(Y, 2, function(y) nnls(P_2, y)$x)
    }
    else if (methods$var == "ols") {
        sigma_var <- solve(t(P_2) %*% P_2) %*% OBOB(t(P_2) %*% Y)
    }
    for (k in 1:(K*(K+1)/2)) {
        diag(Sigma_array[, , k]) <- sigma_var[k, ]
    }
    if (methods$covar == "ols") {
        n <- nrow(P)
        w <- rep(1, n)
    }
    else if (methods$covar == "wls") {
        obs_w <- P_2 %*% sigma_var
    }
    for (i in 1:(p - 1)) {
        for (j in (i + 1):p) {
            if (methods$covar == "wls") {
                w <- sqrt(obs_w[, i] * obs_w[, j])
            }
            P_2_w <- P_2/w
            y_w <- X_centered[, i] * X_centered[, j]/w
            Sigma_array[j, i, ] <- Sigma_array[i, j, ] <- solve(t(P_2_w) %*% 
                P_2_w) %*% t(P_2_w) %*% y_w
        }
    }
    return(lapply(1:(K*(K+1)/2), function(k) Sigma_array[, , k]))
}


make_plots <- function(est_list, 
                       type = 'var',
                       gene_set_names,
                       to_save,
                       var_ls, 
                       method_list=list('CSNet', 'CSNet-cross')){
    require(latex2exp)
    require(dplyr)
    if(!to_save) options(repr.plot.width = 20, repr.plot.height = 6)
    
    if(to_save){
        pdf(sprintf('figures/%s_%s_vs_%s_%s.pdf', gene_set_names, method_list[1], method_list[2], var_ls),
       #width = 10, height = 10)
            width = 5*4, height = 6)
        
    }
    par(mfrow = c(2,2), cex = 1.5)
    notation <- ifelse(type == 'var', "$|\\hat{\\sigma}_{jj}|$", "$ |\\hat{\\sigma}_{jj'}| $")
    df_list <- list()
    for(ct in major_cts){
        if(type == 'var'){
            x <- (est_list[[1]][[ct]] %>% diag) %>% abs %>% log
            y <- (est_list[[2]][[ct]] %>% diag) %>% abs %>% log
        }else{
            uinds <- upper.tri(est_list[[1]][[1]])
            x <- est_list[[1]][[ct]][uinds] %>% abs %>% log
            y <- est_list[[2]][[ct]][uinds] %>% abs %>% log
        }
        #print(summary(x))
        #print(summary(y))
        df_list[[ct]] <- matrix(c(x, y), byrow = F, ncol = 2)
        #plot(x, y,
        # xlab = TeX(sprintf(r'(log %s, %s)', notation, method_list[[1]])), 
        # ylab = TeX(sprintf(r'(log %s, %s)', notation, method_list[[2]])),
        #main = sprintf('%s, cor=%.2f', ct,
        #               cor(x, y, method = 'pearson')))
        #abline(0,1,col = 'red', lty = 'dashed')
    }
    df <- do.call(rbind, df_list) %>% as.data.frame
    colnames(df) <- c('x', 'y')
    ct_cors <- sapply(major_cts, function(ct) cor(df_list[[ct]][,1], df_list[[ct]][,2]))
    ct_titles <- sapply(1:4, function(i) sprintf('%s \ncor=%.2f', major_cts[i], ct_cors[i]))
    df$ct <- rep(ct_titles, each = nrow(df_list[[1]]))
    g <- ggplot(df) +
    geom_point(aes(x = x, y = y), alpha = ifelse(type == 'covar', 0.2, 1)) +
    labs(x = ifelse(gene_set_names == 'AD risk genes', TeX(sprintf('log %s, %s', notation, method_list[[1]])), ''),
         
        y = TeX(sprintf('log %s, %s', notation, method_list[[2]])),
         title = gene_set_names
         #title = NULL # remove title since there's caption
         #title = ifelse(type == 'var', 'Cell-type-specific variances', 'Cell-type-specific covariances')
        ) +
        # title = gsub('_', ' ', gene_set_names)) +
    geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 'dashed') +
    theme_classic(base_size = 30) +
    facet_wrap(~ct, nrow = 1)
    print(g)
    
    if(to_save) dev.off()
}

compare_two_CSNet <- function(gene_set_names, var_ls = 'ols', to_save = F){
    require(stringr)
    gene_set <- load_genes(gene_set_names)
    #tmp <- strsplit(gene_set_names, '_')[[1]]
    #gene_set_names <- paste0(str_to_title(tmp[2]), ' ', tolower(tmp[3]))

    ## 2. construct bulk data
    gt_cor_list <- gt_cov_list <- list()
    pb_bulk_X <- matrix(0, nrow = nrow(main_ct_prop_tab), ncol = length(gene_set))
    for(ct in major_cts){
        # bulk data with weights given by observed cell type proportions
        st_exp <- t(data_list[[ct]]$scaled_data_cons_depth[gene_set, ]) %>% as.matrix
        pb_bulk_X <- pb_bulk_X +  st_exp * P[, ct]
        # oracle / ground truth estimates
        gt_cor_list[[ct]] <- cor(data_list[[ct]]$scaled_data_cons_depth[gene_set, ] %>% t)
        gt_cov_list[[ct]] <- cov(data_list[[ct]]$scaled_data_cons_depth[gene_set, ] %>% t)
    }
    ## 3. obtain estimates from two versions of CSNet
    # fix covariance to be estimated by wls
    # without cross-cell-type
    dCSNet_without_cov <- mom_ls(P, pb_bulk_X, methods = list(var=var_ls, covar = 'wls'))
    dCSNet_without <- lapply(dCSNet_without_cov, function(x) get_cor_from_cov(x))
    names(dCSNet_without) <- names(dCSNet_without_cov) <- major_cts
    # with cross-cell-type
    dCSNet_with_cov_full <- mom_ls_with_cross(P, pb_bulk_X, methods = list(var=var_ls, covar = 'wls'))
    dCSNet_with_cov <- lapply(c(1,5,8,10), function(i) dCSNet_with_cov_full[[i]])
    names(dCSNet_with_cov) <- major_cts
    # save plots
    # (1) compare variances
    if(var_ls == 'ols') make_plots(list(dCSNet_without_cov, dCSNet_with_cov), 'var', 
                                   gene_set_name_titles[gene_set_names], to_save, var_ls)
    if(var_ls == 'nnls') make_plots(list(dCSNet_without_cov, dCSNet_with_cov), 'covar', 
                                    gene_set_name_titles[gene_set_names], to_save, var_ls)
    # return(list(dCSNet_without_cov, dCSNet_with_cov))
}


