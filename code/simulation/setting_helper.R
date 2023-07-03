#
# setting_helper.R
#
# generate parameters used for simulation settings

gen_sim_setting <- function(cor_model,
	n, p, log_var, 
	K, betas, 
	AR_MA_ctrl = list(rhos = c(0.9, 0.9)),
	K_4_ctrl = list(equal_strength = T),
	real_ctrl = list(real_threshold = 'th_0.6'),
	NB_exper_pars = list(Tsize=200, S=6e+07), 
	verbose=F,
	equal_sigma_sq = F){

	if(cor_model != 'dcSBM'){
		sim_setting <- gen_sim_setting_AR_MA(n, p, log_var, 
			cor_model, AR_MA_ctrl$rhos,
			K, betas,
			K_4_ctrl,
			real_ctrl,
			NB_exper_pars, 
			verbose,
			equal_sigma_sq)
	}else{
		
	}
}

# generate simulation parameters for the AR(1)/MA(1) model with exact NB simulation
gen_sim_setting_AR_MA <- function(n, p, log_var, 
	cor_model, rhos,
	K, betas, 
	K_4_ctrl,
	real_ctrl,
	NB_exper_pars = list(Tsize=200, S=6e+07), 
	verbose = F,
	equal_sigma_sq = F){

	source('NB_helper.R')
	print(sprintf('Generate simulation parameters for a K=%i, log_var=%.1f, %s model', 
		K, log_var, cor_model))
	
	# set dimension and ids for cell-type-specific gene sets
	#if(cor_model != 'real_data'){
	if(K == 2 | K == 4){
		cor_p <- round(p / (K + 1))
		sub_cl <- lapply(1:K, function(k) (((k-1)*cor_p+1) :(k*cor_p)))
	}else if(K == 10){
		K0 <- 4
		cor_p <- round(p / (K0 + 1))
		sub_cl <- lapply(1:K0, function(k) (((k-1)*cor_p+1) :(k*cor_p)))
		for(k in 5:8) sub_cl[[k]] <- sub_cl[[k-4]]
		sub_cl[[9]] <- sub_cl[[1]]
		sub_cl[[10]] <- sub_cl[[2]]
	}#}else{
	#	cor_p <- p
	#	sub_cl <- lapply(1:K, function(k) 1:p)
	#}
	# simulate and fix cell type proportions
	set.seed(1)
	if(K == 2){
		beta1 <- betas[1]
		beta2 <- betas[2]
		pi_ct1 <- rbeta(n, beta1, beta2)
		pi_m <- matrix(c(pi_ct1, 1-pi_ct1), ncol = 2)
		colnames(pi_m) <- paste('cell type', c(1,2))
	}else if(K > 2){
		require(rBeta2009)
		pi_m=rdirichlet(n, betas)
		colnames(pi_m) <- paste('cell type', 1:K)
	}

	# set expression variance
	if(K == 2){
		sigma_sq_1 <- rep(exp(log_var), p)
		sigma_sq_2 <- rep(exp(log_var), p)
		#if(cor_model != 'real_data'){
		if(!equal_sigma_sq){
			# to demonstrate the confounding effect of cell type proportions
			# we consider a gene cluster with independent expression yet different mean expression levels
			# in two cell types
    			sigma_sq_2[(2*cor_p+1): p] <- rep(exp(log_var+1), p-2*cor_p) 
			#}else{
			#	sigma_sq_2[(2*round(p/3)+1): p] <- rep(exp(log_var+1), p-2*round(p/3))
			#}
		}
		sigma_sq_star_list <- list(sigma_sq_1, sigma_sq_2)
	}else if(K == 4 | K == 10){
		sigma_sq_star_list <- lapply(1:K, function(k) rep(exp(log_var), p))
    		less_abun_celltypes <- (1:K)[!(1:K) %in% which.max(betas)]
    		if(!K_4_ctrl$equal_strength){
      			sigma_sq_star_list[[4]][sub_cl[[4]]] <- rep(exp(log_var+2), cor_p)
    			#sigma_sq_star_list[[4]] <- rep(exp(log_var+3), p)
		}
    		for(k in less_abun_celltypes){
      			sigma_sq_star_list[[k]][(4*cor_p+1):p] <- rep(exp(log_var+1), p-4*cor_p)
    		}
	}

	# generate correlation matrices
	R_star_list <- list()
	NB_par_list <- Sigma_star_list <- R_star_list <- list()
	for(k in 1:K){
		cor_ind <- sub_cl[[k]]
		NB_par_list[[k]] <- simu_NB_cov(cor_p, rho=rhos[k], model=cor_model,
                                    Tsize=NB_exper_pars[['Tsize']], S=NB_exper_pars[['S']],
                                    use_approx=T,
                                    var_vec = sigma_sq_star_list[[k]][cor_ind],
				    k=k,
				    real_threshold=real_ctrl$real_threshold)
        	Sigma_star <- diag(sigma_sq_star_list[[k]])
        	Sigma_star[cor_ind, cor_ind] <- NB_par_list[[k]][['cov']]
    		Sigma_star_list[[k]] <- Sigma_star
    		R_star_list[[k]] <- Sigma_star / outer(sqrt(sigma_sq_star_list[[k]]),
    			sqrt(sigma_sq_star_list[[k]]))
    }

    # generate cell-type-specific mean expressions
    mu_star_list <- list()
	for(k in 1:K){
    		gene_mu_vec <- numeric(p)
    		for(j in 1:p){
      			phi <- set_phi(sigma_sq_star_list[[k]][j],
         	            Tsize=NB_exper_pars$Tsize, S=NB_exper_pars$S)
      		gene_mu <- NB_exper_pars$S * NB_exper_pars$Tsize * phi
      		gene_mu_vec[j] <- gene_mu
    		}
    		mu_star_list[[k]] <- gene_mu_vec
  	}

  	sim_setting <- list('mu' = mu_star_list,
                       	'sigma_sq' = sigma_sq_star_list,
                       	'R' = R_star_list,
                       	'props' = pi_m,
                       	'NB_par' = NB_par_list,
                       	'cor_p' = cor_p,
                       	'sub_cl' = sub_cl,
                       	'cor_model' = cor_model)
  	return(sim_setting)
}


# generate simulation parameters for the dcSBM model with copula

gen_sim_setting_dcSBM <- function(n, p, log_var, 
                                K, betas, 
                                K_4_ctrl = list(equal_strength = T),
                                NB_exper_pars = list(Tsize=200, S=6e+07)){
	# simualte cell type prorpotions
	require(rBeta2009)
	require(matrixcalc)
	set.seed(1)
	pi_m=rdirichlet(n, betas)
	colnames(pi_m) <- paste('cell type', 1:K)

	# generate correlation matrix from the dcSBM two-block model

}
