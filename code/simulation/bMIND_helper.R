#
# bMIND_helper.R
#
# estimate cell-type-specific co-expressions by 
# the sample correlations of
# cell-type-specific and sample-specific expressions estimated by bMIND


coexp_by_bMIND <- function(data_list, prior_info = 'informative'){
	# install.packages("remotes")
	# remotes::install_github("randel/MIND")
	# remotes::install_github('ChangSuBiostats/bMIND_remove_parallel')
	require(MIND)
	#source('bMIND/MIND/R/bmind_func.r')
	#source('bMIND/MIND/R/bmind_de.r')
	#source('bMIND/MIND/R/deconv_em_func.r')
	p <- ncol(data_list$data$X)
	n <- nrow(data_list$data$X)
	K <- ncol(data_list$data$P)

	# hyper-parameters in bMIND
	if(prior_info == 'informative'){
		# default in their implementation
		exp_prior_var <- 0.5
		nu <- 50
	}else{
		exp_prior_var <- 10000000000
		nu <- 0
	}

	# -
	# format the data as input to bMIND
	# -
	bulk <- t(data_list$data$X)
	rownames(bulk) <- paste0('Gene_', 1:p)
	colnames(bulk) <- paste0('Sample_', 1:n)

	ref <- do.call(cbind, lapply(1:K, function(k) t(data_list$ct_specific_data[[k]])))
	rownames(ref) <- rownames(bulk)
	colnames(ref) <- sapply(1:K, function(k) paste(colnames(bulk), sprintf('ct%i', k), sep = '_'))
	ref_meta <- data.frame(sample = rep(colnames(bulk), K),
		cell_type = rep(paste0('ct_', 1:K), each = n))
	rownames(ref_meta) <- colnames(ref)

	# normalize both bulk and cell-sorted data to be on CPM scale
	bulk <- edgeR::cpm(bulk)
	ref <- edgeR::cpm(ref)

	# -
	# run bMIND
	# -

	# evaluate prior
	prior = get_prior(sc = ref, meta_sc = ref_meta)

	# obtain true cell type proportions
	frac <- data_list$data$P
	rownames(frac) <- colnames(bulk)
	colnames(frac) <- paste0('ct_', 1:K)

	# evaluate bMIND estimates
	deconv = bMIND(bulk = log2(bulk+1),
                   frac = frac,
                   sample_id = colnames(bulk),
                   ncore = NULL,
                   profile = prior$profile,
                   covariance = prior$covariance,
                   nu = nu,
                   V_fe = diag(exp_prior_var, ncol(prior$profile)))
	R_est <- lapply(1:K, function(k){
				tmp <- cor(t(deconv$A[,k,]))
				if(anyNA(tmp)){
					# post-process the matrix if estimated cell-type-specific expression levels are constant across samples
					n_cons <- apply(tmp, 1, function(x) all(is.na(x))) %>% length
					print(sprintf('%i genes were estimated as constant by bMIND in cell type %i', n_cons, k))
					tmp[is.na(tmp)] <- 0
					diag(tmp) <- 1
		   }
		   return(tmp)
		   })
	return(R_est)
}

