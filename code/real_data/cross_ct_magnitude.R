# organized from check_cross_ct_indep_ass_400.ipynb

# -
# load packages and codes
# -

library(Matrix)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

# load helper functions for this analysis
source('cross_ct_helper_functions.R')

# load pseudo-bulk data
data_dir <- '/gpfs/gibbs/pi/zhao/cs2629/ROSMAP/GeneExpression/snRNAseq_DLPFC_experiment2/pseudo_bulk'
# sample size per cell type
annot <- read.csv('/gpfs/gibbs/pi/zhao/cs2629/ROSMAP/GeneExpression/snRNAseq_DLPFC_experiment2/cell-annotation.csv')

ct_props <- table(annot$cell.type) / nrow(annot)
major_cts <- c('Excitatory Neurons', 'Oligodendrocytes',
              'Astrocyte', 'Microglia')

data_list <- lapply(major_cts, function(ct) load_data(ct))
names(data_list) <- major_cts

# gene sets of interest
AD_genes <- load_genes('AD_risk_genes')
Ex_genes <- load_genes('GOCC_EXCITATORY_SYNAPSE')
Oli_genes <- load_genes('GOCC_MYELIN_SHEATH')
Ast_genes <- load_genes('GOBP_ASTROCYTE_DIFFERENTIATION')

# full gene names
features <- read.table(sprintf('%s/%s_features.txt', data_dir, major_cts[1]))[[1]]

# select genes that have expression levels higher than rank 10000
# in all cell types
top_inds <- list()
for(top_n in c(10000)){
    top_k_inds <- sapply(major_cts, function(ct){
        data_list[[ct]]$gene_summary$mean_exp_order < top_n
    })
    top_inds[[as.character(top_n)]] <- apply(top_k_inds, 1, function(x) all(x))
}

# -
# plot the distribution of log(|within|/|between|)
# for four gene sets
# -

plot_within_versus_between(AD_genes, T, 'AD_risk_genes')
plot_within_versus_between(Ex_genes, T, 'Excitatory synapse')
plot_within_versus_between(Oli_genes, T, 'Myelin sheath')
plot_within_versus_between(Ast_genes, T, 'Astrocyte differentiation')
