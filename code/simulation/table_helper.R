#
# table_helper.R
#
# helper functions for printing tables
print_a_setting <- function(error_list, coexp_methods, n, p, metric_index = 1:4){
  n_methods <- length(coexp_methods)
  if(p <= 100){
    cat('\\\\\\hline') 
    cat('\n')
    cat(sprintf('\\multirow{%i}{*}{%i} & \\multirow{%i}{*}{\\begin{tabular}[c]{@{}l@{}}$%i$\\end{tabular}}',
        4 * n_methods, n, 2 * n_methods, p))
  }else{
    cat('\\\\\\cline{2-11}')
    cat('\n')
    cat(sprintf('& \\multirow{%i}{*}{\\begin{tabular}[c]{@{}l@{}}$%i$\\end{tabular}}',
        2 * n_methods, p))
  }
  cat('\n')
  t <- 1
  for(coexp_method in coexp_methods){
    print_a_method(error_list[[coexp_method]], coexp_method, t == 1, metric_index)
    t <- t + 1
  }
  # cat('\\\\\\hline')
}

print_a_sensitivity_setting <- function(error_list, coexp_methods, kappa, b, rho_cor, rMSE){
  cat('\n')
  if(b == -0.4){
    cat('\\\\\\hline')
    cat(sprintf('\\multirow{%i}{*}{%.1f}',
        2 * 5, kappa)) #, b, rho_cor, rMSE))
  }else{
    cat('\\\\\\cline{2-12}')
    #cat(sprintf('& %.1f & %.2f & %.2f',
    #    b, rho_cor, rMSE))
  }
  cat('\n')
  cat('\n')
  for(coexp_method in coexp_methods){
    print_a_sensitivity_method(error_list[[coexp_method]], b, rho_cor, rMSE)
  }
}

print_a_method <- function(error_mat, coexp_method, first_method, metric_index = 1:4){
  K <- length(error_mat)
  if(K == 10){
    # if K=10, visualize only the first 4 cell types
    K <- 4
  }
  if(!first_method){
    cat('\\\\\\cline{3-11}')
    cat('\n')
    cat(sprintf('&& \\texttt{%s}', coexp_method))
  }else{
    cat(sprintf('& \\texttt{%s}', coexp_method))
  }
  cat('\n')
  #for(j in 1:ncol(error_mat[[1]])){
  for(k in 1:K){
    for(j in metric_index){
    #for(j in 1:ncol(error_mat[[1]])){
    err_mean <- mean(error_mat[[k]][,j])
    err_sd <- sd(error_mat[[k]][,j])
    if(coexp_method == 'CSNet'){
      cat('& \\begin{tabular}[c]{@{}c@{}}',
          ifelse(is.na(err_mean), '{-}', sprintf('\\bf{%.2f}', err_mean)),
          '\\\\\\small{',
          ifelse(is.na(err_sd), '{-}', sprintf('\\bf{(%.2f)}', err_sd)),
          '}\\end{tabular}', sep='')
    }else{
      cat('& \\begin{tabular}[c]{@{}c@{}}',
          ifelse(is.na(err_mean), '{-}', sprintf('%.2f', err_mean)),
          '\\\\\\small{',
          ifelse(is.na(err_sd), '{-}', sprintf('(%.2f)', err_sd)),
          '}\\end{tabular}', sep='')
    }
    }
  }
  cat('\n')
}

print_a_sensitivity_method <- function(error_mat, b, rho_cor, rMSE){
  K <- length(error_mat)
  cat('\n')
  cat(sprintf('& %.1f & %.2f & %.2f', b, rho_cor, rMSE))
  cat('\n')
  #for(j in 1:ncol(error_mat[[1]])){
  for(k in 1:K){
    for(j in 1:ncol(error_mat[[1]])){
    err_mean <- mean(error_mat[[k]][,j])
    err_sd <- sd(error_mat[[k]][,j])
      cat('& \\begin{tabular}[c]{@{}c@{}}',
          ifelse(is.na(err_mean), '{-}', sprintf('%.2f', err_mean)),
          '\\\\\\small{',
          ifelse(is.na(err_sd), '{-}', sprintf('(%.2f)', err_sd)),
          '}\\end{tabular}', sep='')
    }
  }
  cat('\n')
}
