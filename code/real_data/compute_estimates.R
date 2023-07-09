
library("optparse")
option_list = list(
  make_option(c("--i_geneset"), type="integer", default=1,
              help="index of GO gene set", metavar="integer")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

i_geneset <- opt$i_geneset

################
# Load packages, codes and data
###############

# -
# packages
# -
library(dplyr)
library(nnls)

# for ordering genes in the heatmap
library(WGCNA)

# for handling single cell data
library(Seurat)

# for visualization
library(ggplot2)
library(grid)
library(gridExtra)

# load bMIND and ENIGMA
library(MIND)

# -
# codes
# -
source('/home/cs2629/project/CSNet/JASA_codes/CSNet_analysis/code/mom_ls.R')
source('/home/cs2629/project/CSNet/JASA_codes/CSNet_analysis/code/helper_function.R')
source('/home/cs2629/project/CSNet/JASA_codes/CSNet_analysis/code/cross_validation.R')
source('/home/cs2629/project/CSNet/JASA_codes/CSNet_analysis/code/simulation/visualization_helper.R')
source('real_data_helper_functions.R')

# -
# parameters
# -

# cell types of interests
cts <- c('Ex', 'Oli', 'Ast', 'Mic')
ct_full_names <- c('Excitatory neuron', 'Oligodendrocyte', 'Astrocyte', 'Microglia')

# gene set of interests
genesets <- c('GOCC_EXCITATORY_SYNAPSE', 
	      'GOCC_MYELIN_SHEATH', 
	      'GOBP_ASTROCYTE_DIFFERENTIATION',
	      'AD_risk_genes')

# dimension for bulk heatmaps
bulk_hm_dim_list <- list(c(1.3, 1.2), c(1.3, 1.2), c(1.6, 1.7), c(1.4, 1.6)) # c(1.4, 1.72)
names(bulk_hm_dim_list) <- genesets
plot_width <- 4


# -
# directories
# -

# This folder contains RNA-seq data from the ROSMAP project, which are under controlled access.
# Application for the access can be found here:
# https://www.synapse.org/#!Synapse:syn3388564
rosmap_data_dir <- '/gpfs/gibbs/project/fan_zhou/cs2629/CSNet/real_data/covest-real-data/ROSMAP'

# This folder contains single nucleus data from the ROSMAP project, which are under controlled access.
# Application for the access can be found here:
#  https://www.synapse.org/#!Synapse:syn21261143
rosmap_sc_data_dir <- '/gpfs/gibbs/pi/zhao/cs2629/ROSMAP/GeneExpression/snRNAseqPFC_BA10'

# This folder stores intermediate data results that cannot be shared
# due to the same reason as above.
output_dir <- 'output/'
figure_dir <- 'figures/coexp_estimates'

# Data in this folder are available with the github repo
saved_data_dir <- '../../data/ROSMAP'


# -
# load preprocessed bulk data for the specific gene set
# -
geneset <- genesets[i_geneset]
print(sprintf('Compute co-expression estimates for %s gene set', geneset))

data_res <- readRDS(sprintf('%s/%s_data.rds', output_dir, geneset))
full_data_list <- data_res$full_data_list
data_list <- full_data_list[[1]]
genes <- full_data_list[[2]]
bulk_cl_hc <- data_res$bulk_cl_hc

all_cts <- colnames(data_list$P)

###############
# compute estimates
###############

# -
# CSNet estimates
# -

# NNLS for variances and WLS for covariances
methods_settings <- list(var='nnls', covar='wls')
# threshold for small weights in WLS
th <- 0.01

dCSNet_cov_est <- mom_ls(data_list$P, data_list$X, methods_settings, th)
dCSNet_est <- lapply(dCSNet_cov_est, function(x) get_cor_from_cov(x))


# d-CSNet
plot_four_heatmaps(dCSNet_est, bulk_cl_hc, geneset, 'd-CSNet', cw = NULL)
                                    
# CSNet

## regular cross validation
cv_res <- cross_validation(data_list,
                           coexp_method = 'CSNet',
                           list(methods_settings = methods_settings, th = th),
                           ct_names = colnames(data_list$P),
                           to_plot = F,
                           seed = 1,
                          n_splits = 1) 

CSNet_est <- lapply(1:length(dCSNet_est), function(i){
                            generalized_th(dCSNet_est[[i]], cv_res$th[i], 'SCAD', F) })
plot_four_heatmaps(CSNet_est, bulk_cl_hc, geneset, 'CSNet', cw = NULL)
saveRDS(CSNet_est, sprintf('%s/%s_CSNet_est.rds', saved_data_dir, geneset))

if(i_geneset == 4){
	## 1 se
	cv_res <- cross_validation(data_list,
				   coexp_method = 'CSNet',
				   list(methods_settings = methods_settings, th = th),
				   ct_names = colnames(data_list$P),
				   to_plot = F,
				   seed = 1,
          n_splits = 10)
	print(cv_res$th)
	print(cv_res$th_1se)
	CSNet_est <- lapply(1:length(dCSNet_est), function(i){
				    # apply the 1-se rule only to the less abundant cell types,
				    # including oligodendrocyte, astrocyte and microglia (and other rare cell tyeps)
				    th <- ifelse(i == 1, cv_res$th[i], cv_res$th_1se[i])
				    print(th)
				    generalized_th(dCSNet_est[[i]], th, 'SCAD', F) 
	  })
	plot_four_heatmaps(CSNet_est, bulk_cl_hc, geneset, 'CSNet_1se_for_less_abundant', cw = NULL)
	saveRDS(CSNet_est, sprintf('%s/%s_CSNet_est_1se_for_less_abundant.rds', saved_data_dir, geneset))
}

# -
# bulk estimates
# -
bulk_coexp <- cor(data_list$X)

g <- plot_heatmap(bulk_coexp[bulk_cl_hc$reordering, bulk_cl_hc$reordering],
             bulk_cl_hc$cluster,
             'Bulk brain',
             legend = T,
             annotation_legend = T,
             annotation_names = F,
             cw = 4,
             silent = T)[[4]]
g <- grid.arrange(grobs = list(g), nrow = 1)
ggsave(sprintf('%s/%s_bulk_sample_cor.pdf', figure_dir, geneset), g,
       width = plot_width * bulk_hm_dim_list[[geneset]][1],
       height = plot_width * bulk_hm_dim_list[[geneset]][2])


# -
# s-bMIND estimates
# -

# load bMIND priors evaluated in compute_bMIND_priors.R
prior <- readRDS('output/bMIND_prior.rds')

# split genes based on whether or not it have valid priors
# inferred from bMIND
ind_in_bMIND_bulk <- which(genes %in% rownames(prior$profile))
ind_in_bMIND <- match(genes[ind_in_bMIND_bulk], rownames(prior$profile))
print(sprintf('%i genes have valid priors inferred from bMIND',
              length(ind_in_bMIND_bulk)))

# obtain bMIND estimate (dense)
bMIND_est <- real_data_run_bMIND(data_list, ind_in_bMIND_bulk, ind_in_bMIND)


# use cross validation to evaluate s-bMIND
tuning_result <- cross_validation(data_list,
                                 coexp_method = 'bMIND',
                                 method_pars = list(ind_bulk = ind_in_bMIND_bulk,
                                                    ind_profile = ind_in_bMIND),
                                 to_plot = F,
				 seed = 1,
				 n_splits = 1)

# obtain s-bMIND estimates
# sbMIND_est_1se
sbMIND_est_min_th <- list()
for(k in 1:length(all_cts)){
	# not consider NA entries when doing thresholding
    tmp <- bMIND_est[[k]]
    na_inds <- is.na(tmp)
    tmp[na_inds] <- 0
    diag(tmp) <- 1

    tmp_min_th <- generalized_th(tmp, tuning_result$th[k], 'SCAD', F)

    tmp_min_th[na_inds] <- NA
    diag(tmp_min_th) <- 1
    sbMIND_est_min_th[[k]] <- tmp_min_th
}

plot_four_heatmaps(bMIND_est, bulk_cl_hc, geneset, 'bMIND', cw = NULL)
plot_four_heatmaps(sbMIND_est_min_th, bulk_cl_hc, geneset, 's-bMIND', cw = NULL)


# -
# single cell based
# -
sc_est <- list()
for(ct in cts){
    pseudo_bulk_data <- load_bulk_sc_data(genes, ct)
    # compute sample correlations
    sc_sample_cor <- cor(pseudo_bulk_data %>% t)
    # use cross validation to select thresholding parameters
    tuning_result <- cross_validation(list(X = pseudo_bulk_data %>% t))
    # generalized thresholding - SCAD
    tmp <- sc_sample_cor
    na_ind <- is.na(tmp)
    tmp[na_ind] <- 0
    diag(tmp) <- 1

    sc_est[[ct]] <- generalized_th(tmp, tuning_result$th, 'SCAD', F)

    sc_est[[ct]][na_ind] <- NA
    diag(sc_est[[ct]]) <- 1
}

plot_four_heatmaps(sc_est, bulk_cl_hc, geneset, 'single_cell', cw = NULL)



