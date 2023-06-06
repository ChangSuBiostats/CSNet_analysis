# simulate multivariate NB random variables
# using a 3-step data generating process model
# 
# Chang Su, c.su@yale.edu

# ----------------
# helper functions 
# ----------------

# derive ground truth variances based on parameters
eval_var <- function(S, Tsize, phi){
  return(c(pred = S^2 * Tsize * phi^2 + S * Tsize * phi,
           pred_approx = S^2 * Tsize * phi^2))
}

eval_cor <- function(S, Tsize, phi, T12){
  return(c(pred = T12 / sqrt(((1 / (S * phi[1]) + 1) * Tsize) * 
                               ((1 / (S * phi[2]) + 1) * Tsize)),
           pred_approx = T12 / Tsize))
}

# ------------------------------------------------------
# set NB parameters using var, cov and cor specification
# ------------------------------------------------------

# set parameters based on var/cov
set_phi <- function(sigma_sq, Tsize=200, S=6e+07, use_approx=F){
  if (use_approx){
    phi <- 1/S * sqrt(sigma_sq / Tsize)
  }else{
    phi <- 1/S * sqrt(sigma_sq / Tsize + 1/4) - 1/(2*S)
  }
  return(phi)
}

# set T12 based on correlation
# caveat: this function implicitly assume T to be constant
#         across genes
set_T_shared <- function(rho, phi, Tsize, S, use_approx=F){
  if (use_approx){
    T_shared <- rho * Tsize
  }else{
    T_shared <- rho * Tsize * sqrt((1/(S * phi[1]) + 1) * (1/(S * phi[2]) + 1))
  }
  return(round(T_shared))
}

# get |A_i \cap A_j| for any (i,j)
# for i,j=1,\dots,p
get_T_shared <- function(cor_m, phi_vec, S, Tsize=200,use_approx=F){
  # if use_approx = T
  # use approximation and let Tsize proportional to rho
  # with the drawback that, the generated covariance matrix
  # is not exactly AR(1) when converting to correlations for NB
  # but it is for Gammas
  if (use_approx){
    size_m <- round(Tsize * cor_m)
  }else{
    adj <- outer(1/(S * phi_vec) + 1, 1/(S * phi_vec) + 1)
    # rounding to get integers for #components to share
    size_m <- round(cor_m * Tsize * sqrt(adj))
    diag(size_m) <- Tsize
  }
  return(size_m)
}

# ---------------------------------------
# simulate multivariate NB 
# by sharing exponential 'blocks'
# ---------------------------------------

map_to_vec <- function(i, j, p){
  return(i + (j-1) * (j-2)/2)
  # return((i - 1) * (p - i /2) + (j - i))
  # return(i + (j - 1) * (j - 2) / 2)
}

# step 1 and 2:
# generate exponential building blocks and
# aggregate over overlapped and indpt exponential blocks to create
# correlated Gamma r.v.
gen_cor_gammas <- function(size_m, phi_vec, Tsize=200){
  m <- ncol(size_m)
  size_inds <- cumsum(size_m[upper.tri(size_m)])
  # simulate shared blocks
  blocks <- rexp(max(size_inds), rate=1)
  # simulate multivariate Gamma r.v. p
  p <- numeric(m)
  for(i in 1:m){
    # get flattened indexes for the i-th r.v.
    inds <- numeric(m - 1)
    t <- 1
    for(j in 1:m){
      if(j < i){
        inds[t] <- map_to_vec(j, i, m)
        t <- t + 1
      }else if(j > i){
        inds[t] <- map_to_vec(i, j, m)
        t <- t + 1
      }
    }
    # print(inds)
    # get actual indexes for the blocks
    b_inds <- sapply(inds, function(ind){
      if(ind > 1){
        if(size_inds[ind-1] < size_inds[ind]){
          (size_inds[ind-1] + 1):size_inds[ind]
        }
      }else{
        if(size_inds[ind] > 0){
          1:size_inds[ind]
        }
      }
    })
    # print(b_inds)
    b_inds <- unlist(b_inds)
    p[i] <- sum(blocks[b_inds]) * phi_vec[i] + 
      sum(rexp(as.integer(Tsize-length(b_inds)), rate=1)) * phi_vec[i]
  }
  return(p)
}

# or this function:
# specialized for approxiamting AR1 model
gen_cor_gammas_AR1 <- function(size_m, phi_vec, Tsize=200){
  require(dplyr)
  
  ngenes <- nrow(size_m)
  stop_ind <- max(which(size_m[1,] >= 1)) - 1
  
  # make a size matrix
  # to keep track of the accumulation of correlations
  # note:
  # the 'blocks' are accumlated from the right to the left horizontally
  # while they are independent from the top to the bottom
  # i.e. A_{i,j+1} \subset A_{i,j}
  #      A_{i,j} \cap A_{i+1,j} = {} (empty set)
  cumu_size_m <- matrix(0, ngenes, ngenes)
  for(i in 1:(ngenes - 1)){
    for(j in 1:stop_ind){
      col_ind <- min(i + j, ngenes)
      if(i > 1){
        cumu_size_m[i, col_ind] <- size_m[i, col_ind] - size_m[i-1, col_ind]
      }else{
        cumu_size_m[i, col_ind] <- size_m[i, col_ind]
      }
    }
  }
  diag(cumu_size_m) <- 0
  # make a vector 
  # to store total number of blocks per row
  row_total <- apply(cumu_size_m, 1, function(c_row) max(c_row)) 
  
  # generate 'building blocks' by row
  # and aggregate to get correlated Gamma random variables
  gamma_obs <- gamma_obs_1 <- gamma_obs_2 <- numeric(length = ngenes)
  
  for(i in 1:(ngenes)){
    if(row_total[i] > 0){
      blocks <- rexp(row_total[i], rate=1)
    }else{
      blocks <- 0
    }
    # for each row, the corresponding Gamma r.v. would 
    # aggregate over all correlated components
    gamma_obs_1[i] <- sum(blocks)
    
    if(i < ngenes){
      for(j in 1:stop_ind){
        # for each column, the Gamma r.v. would
        # aggregate over corresponding overlaped components
        col_ind <- i + j
        if(col_ind <= ngenes){
          if(cumu_size_m[i, col_ind]>=1){
            gamma_obs_2[col_ind] = gamma_obs_2[col_ind] + 
              sum(blocks[1:cumu_size_m[i, col_ind]])
          }
        }
      }
    }
  }
  for(i in 1:ngenes){
    gamma_obs[i] <- gamma_obs_1[i] + gamma_obs_2[i]
  }
  
  # count total number of correlated component per r.v.
  cor_total <- colSums(cumu_size_m) + 
    c(sapply(1:(ngenes-1), function(i) cumu_size_m[i, i + 1]), 0)
  
  #print(row_total[40])
  #print(gamma_obs_1[40])
  #print((cor_total - row_total)[40])
  #print(mean(gamma_obs_2))
  #print(mean(gamma_obs))
  
  
  # check if indeed this AR(1) model can be generated by our model
  if(any(cor_total > Tsize)){
    warning('sum of correlated building blocks exceed the total number per r.v.')
  }
  
  # fill in the 'rest' of each Gamma random variables by 
  # indpt exponential blocks
  for(i in 1:ngenes){
    if(Tsize < cor_total[i]){
      print(i)
      print(cumu_size_m[i,])
      print(cor_total[i])
      print(size_m[i,])
    }
    gamma_obs[i] <- gamma_obs[i] + rexp(Tsize - cor_total[i], rate = 1) %>% sum
  }
  
  #print(mean(gamma_obs))
  # adjust each NB obervation by its phi
  gamma_obs <- gamma_obs * phi_vec
  return(gamma_obs)
}

# step 3:
# simulate NB by Poisson sampling from correlated Gamma r.v.s
gen_cor_NB <- function(gamma_rv, S){
  return(sapply(gamma_rv, function(g) {
    if(g>0){
      rpois(1, S*g)
    }else{
      0
    }
  }))
}

## Alternatively,
## Instead of exact simulation,
## use copula to simulate correlated NB random variables

gen_cor_NB_copula <- function(n, cor_mat, sigma_sq, mu){
        p <- nrow(cor_mat)
        cor_mat_up <- chol(cor_mat)
    copula <- matrix(rnorm(p*n),nrow = p)
    copula <- pnorm(t(cor_mat_up) %*% copula)

    exp_matrix <- matrix(NA, nrow=p, ncol=n)
    rownames(exp_matrix) <- paste0("Gene", 1:p)

    size_vec <- mu^2 / (sigma_sq - mu)
    for (i in 1:p){
      exp_matrix[i,] <- qnbinom(copula[i,], mu = mu[i], size = size_vec[i])
    }
    return(exp_matrix %>% t)
}

# --------------------------------------
# Simulate simulation parameters
# ---------------------------------------

# simulate a sparse cor matrix
# using two steps:
# 1. AR(1) model with specified rho
# 2. Random permutation to remove ordering
gen_cor <- function(p, rho = 0.3, to_permu=F, model='MA', k=1, real_threshold='th_0.6'){
  cor_m <- matrix(0, nrow = p, ncol = p)
  if(grepl('AR', model)){
    # 1. generate a matrix using AR(1) model
    cors <- rho^(0:(p-1))
    if(model == 'AR'){
      for(j in 1:p){
        cor_m[j, j:p] <- cors[1:(p-j+1)]
      }
    }else{
      m <- as.integer(strsplit(model, '_')[[1]][2])
      for(j in 1:p){
        max_ind <- min(p, j+m)
        cor_m[j, j:max_ind] <- cors[1:(max_ind-j+1)]
      }
    }
  }else if(grepl('MA', model)){
    if('MA' == model){
      m <- 1
    }else{
      m <- as.integer(strsplit(model, '_')[[1]][2])
    }
    # 1. generate a matrix using MA(1) model
    for(i in 1:(p-1)){
      cor_m[i, (i+1):min(i+m,p)] <- rho
    }
  }else if(model == 'real_data'){
    if(!grepl('same_struct', real_threshold)){
      ct <- ifelse(k == 1, 'Excitatory Neurons', 'Oligodendrocytes')
      fname <- sprintf('data/%s_p_%i_%s_spcones_cor.txt', ct, p, real_threshold) 
      cor_m <- read.table(fname) %>% as.matrix()
    }else{
      ct <- 'Excitatory Neurons'
      real_threshold_short <- gsub('_same_struct', '', real_threshold)
      fname <- sprintf('data/%s_p_%i_%s_spcones_cor.txt', ct, p, real_threshold_short)
      cor_m <- read.table(fname) %>% as.matrix()
      set.seed(1)
      if(k == 2){
        permu_inds <- sample(p, p, replace = F)
        cor_m <- cor_m[permu_inds, permu_inds]
      }
    }
  }
  # convert into a symmetric matrix with 1 on the diagonal
  # cor_m[upper.tri(cor_m)] <- cors
  cor_m[lower.tri(cor_m)] <- t(cor_m)[lower.tri(cor_m)]
  diag(cor_m) <- 1
  
  # 2. apply random permutation
  if(to_permu){
    permu_ind <- sample.int(p)
    cor_m <- cor_m[permu_ind, permu_ind]
  }
  return(cor_m)
}

vis_cor <- function(cor_m, title = 'Co-expression'){
  heatmap(cor_m, symm = T, Rowv = NA, Colv = NA,
          col = hcl.colors(12, "Reds", rev = TRUE),
          labRow = FALSE, labCol = FALSE,
          main = title)
}

# simulate variances according to a uniform distribution
sim_vars <- function(p, lower, upper, seed){
  set.seed(seed)
  vars <- runif(p, lower, upper)
  return(vars)
}


# --------------------------------------
# Wrapper around all helper functions
# to give a unified entry to simulate NB r.v.
# ---------------------------------------

simu_NB_cov <- function(ngenes, rho, permu_ind, model='AR', to_permu = F, 
                        Tsize=10, S=6e+07, use_approx=T,v_lower=20000, v_upper=40000, seed=1,
                        var_vec=NULL, k=1, real_threshold='th_0.6'){
  print(sprintf('Simulating %s correlation structure under multivariate NB model', model))
  if(is.null(var_vec)){
    var_vec <- sim_vars(ngenes, v_lower, v_upper, seed)
  }
  phi_vec <- set_phi(var_vec, Tsize, S)
  # generate AR(1) correlation pattern, with possible random permutations
  cor_m <- gen_cor(ngenes, rho=rho, to_permu=to_permu, model=model, 
		   k=k, real_threshold=real_threshold)
  #use_approx <- ifelse(model == 'AR', T, F)
  size_m <- get_T_shared(cor_m, phi_vec, S, Tsize, use_approx)
  # true AR(1) covaraince matrix
  cov_m_target <- cor_m * outer(sqrt(var_vec), sqrt(var_vec))
  if(model %in% c('AR', 'MA')){
    # approximated AR(1) covariance matrix
    cov_m_truth <- size_m * outer(phi_vec, phi_vec) * S^2
    diag(cov_m_truth) <- var_vec
    # if use_approx=F, then these two matrices should be the same
  }else{
    cov_m_truth <- cov_m_target
  }
  # apply permutation after the sharing is simulated
  if(to_permu){
    cor_m <- cor_m[permu_ind, permu_ind]
    cov_m_target <- cov_m_target[permu_ind, permu_ind]
    cov_m_truth <- cov_m_truth[permu_ind, permu_ind]
  }
  
  pars <- list(size_m = size_m, rho = rho, phi_vec = phi_vec, Tsize = Tsize, S = S)
  return(list(pars = pars, cov = cov_m_truth, cov_target = cov_m_target, cor = cor_m))
}

# gen_cor_gammas:
# Algorithm 1: sharing is only allowed between two r.v.s
# adv: precise sharing enumeration; use exact formula for evaluating correlations;
#     can simulate any structure 
# disavd: can only accomodate sum of |cor| \leq 1

# gen_cor_gammas_AR1
# Algorithm 2: recursive sharing is allowed to simulate approximate AR(1)
# avd: can approximate AR(1) with any rho
# disadv: only approximately; AR(1) precisly holds for underlying Gammas

simu_NB_obs <- function(n, pars_list, permu_ind, model='AR', seed=1){
  set.seed(seed)
  pars = pars_list$pars
  nb_obs <- gamma_obs <-matrix(0, nrow = n, ncol = ncol(pars$size_m))
  if(model %in% c('MA', 'AR')){
    for(it in 1:n){
      if(model == 'MA'){
        gammas <- gen_cor_gammas(pars$size_m, pars$phi_vec, pars$Tsize)
      }else if(model == 'AR'){
        gammas <- gen_cor_gammas_AR1(pars$size_m, pars$phi_vec, pars$Tsize)
      }
      nb_obs[it,] <- gen_cor_NB(gammas, pars$S)
      gamma_obs[it,] <- gammas
    }
    nb_obs <- nb_obs[, permu_ind]
  }else{
    mu_vec = pars$Tsize * pars$S * pars$phi_vec
    nb_obs <- gen_cor_NB_copula(n, pars_list$cor, diag(pars_list$cov), mu_vec)
  }
  return(nb_obs)
}
