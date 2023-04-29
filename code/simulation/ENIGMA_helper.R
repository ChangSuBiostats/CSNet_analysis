# 
# ENIGMA_helper.R
#
# estimate cell-type-specific co-expressions by
# the sample correlations of
# cell-type-specific and sample-specific expressions estimated by ENIGMA

# Install ENIGMA following https://github.com/WWXkenmo/ENIGMA
#
# install.packages(c("Matrix","S4Vectors","corpcor","MASS","e1071","ggplot2","cowplot","magrittr","purrr","tibble","nnls","doParallel","tidyr","plyr","vctrs","matrixStats"))
# BiocManager::install(c("SingleCellExperiment","scater","Biobase","SummarizedExperiment","sva","preprocessCore"))
# devtools::install_github("WWXKenmo/ENIGMA_test")

coexp_by_ENIGMA <- function(data_list, norm = 'L2', seed = 1){
  require(ENIGMA)
  p <- ncol(data_list$data$X)
  n <- nrow(data_list$data$X)
  K <- ncol(data_list$data$P)

  # -
  # format the data as input to ENIGMA
  # -
  mixture_df <- data_list$data$X
  colnames(mixture_df) <- paste0('Gene_' , 1:p)
  rownames(mixture_df) <- paste0('Sample_', 1:n)
  
  ## Build a reference matrix by taking the cell-type-specific means for each gene
  ## as suggested by 3 in https://github.com/WWXkenmo/ENIGMA/blob/master/vignettes/A-simple-guide-of-ENIGMA.pdf
  aggre_profile <- sapply(data_list$ct_specific_data, function(x) colMeans(x))
  colnames(aggre_profile) <- paste0('ct_', 1:2)
  rownames(aggre_profile) <- colnames(mixture_df)   

  ## use count data to construct the ENIGMA object as suggested by
  ## https://github.com/WWXkenmo/ENIGMA/blob/master/vignettes/A-simple-guide-of-ENIGMA.pdf
  ## both bulk and reference data are on the count scale (consistent with each other)
  egm = create_ENIGMA(bulk = t(mixture_df),
        ref = aggre_profile,
            ref_type = "aggre")
  
  ## Use true cell type proportions for estimating CSE
  # as what we did for CSNet and bMIND
  # egm = get_cell_proportion(egm, method = "RLR")
  frac <- data_list$data$P
  rownames(frac) <- rownames(mixture_df)
  colnames(frac) <- colnames(aggre_profile)
  egm@result_cell_proportion <- frac

  # -
  # run ENIGMA
  # -
  
  ## Extract CSE
  ## Fit two models with two penalty norm on CSE
  set.seed(seed)
  if(norm == 'L2'){
    egm <- ENIGMA_L2_max_norm(egm)
  }else if (norm == 'trace'){
    egm = ENIGMA_trace_norm(egm)
  }

  # As documented in
  # https://htmlpreview.github.io/?https://github.com/WWXkenmo/ENIGMA/blob/master/vignettes/brain_tutorial.html
  # 'ENIGMA would construct CSE into raw gene expression space'
  # We use the CSE in the raw gene expression space to estimate co-expression
  cse <- sce2array(egm, norm_output = FALSE)
  
  R_est <- lapply(1:K, function(k) cor(t(cse[,,k])))
  return(R_est)
}
