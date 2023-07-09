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
# packages and codes
# -
library(dplyr)
library(WGCNA)
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

# power parameter for WGCNA based on manual examination of the scree plot
geneset_powers <- c(20, 8, 7, 4) 
names(geneset_powers) <- genesets

# threshold for small weights in WLS
th <- 0.01

# -
# data
# -

# This folder contains RNA-seq data from the ROSMAP project, which are under controlled access.
# Application for the access can be found here:
# https://www.synapse.org/#!Synapse:syn3388564
rosmap_data_dir <- '/gpfs/gibbs/project/fan_zhou/cs2629/CSNet/real_data/covest-real-data/ROSMAP'

# This folder stores intermediate data results that cannot be shared
# due to the same reason as above.
output_dir <- 'output/'

# Data in this folder are available with the github repo
saved_data_dir <- '../../data/ROSMAP'

# load bulk data
ROSMAP_bulk <- readRDS(sprintf('%s/output/fpkm_unadj_by_batch.rds', rosmap_data_dir))
bulk_genes <- readRDS(sprintf('%s/output/ROSMAP_annotLookup_new.rds', rosmap_data_dir))
# int_meta <- readRDS('output/ROSMAP_interesting_covar.rds')

# load cell type proportions
train_props_df <- read.table(sprintf('%s/CIBERSORTx/output/train/res_ROSMAP_sc_ref_all_group_5000//CIBERSORTx_Adjusted.txt',
                                     rosmap_data_dir), header = T)
train_props <- train_props_df[,2:9] %>% as.matrix()
dim(train_props)

all_cts <- colnames(train_props)[-which(colSums(train_props) == 0)]

###############
# load data for the specific gene set & obtain gene ordering via WGCNA
###############

if(i_geneset %in% 1:3){
	# -
	# load data
	# -
	geneset <- genesets[i_geneset]
	print(sprintf('Generate data for %s gene set', geneset))
	full_data_list <- load_GO_geneset(geneset, T)
	data_list <- full_data_list[[1]]
	genes <- full_data_list[[2]]
	# -
	# obtain gene clustering results via WGCNA
	# -
	bulk_cl_hc <- run_wgcna(data_list, geneset, geneset_powers[geneset])
}else{
	# load AD genes generated via obtain_AD_gene_list.R
	avail_AD_genes <- read.table(sprintf('%s/avail_AD_genes.txt', saved_data_dir))[[1]]
	geneset <- 'AD_risk_genes'
	data_list_raw <- load_data(avail_AD_genes)
	full_data_list_raw <- list(data_list_raw, avail_AD_genes)
	# in the example of AD risk genes,
	# for better interpretation of the results,
	# we removed 3 genes that cannot be assigned to any clusters in WGCNA.
	bulk_cl_hc_raw <- run_wgcna(data_list_raw, geneset, geneset_powers[geneset], to_plot=F)
	# remove genes that were not assigned to any clusters (cluster 0)
	print(sprintf('Number of unassigned genes: %i', length(bulk_cl_hc_raw$cluster[[1]])))
	genes <- avail_AD_genes[-bulk_cl_hc_raw$reordering[1:length(bulk_cl_hc_raw$cluster[[1]])]]

	new_cluster <- list()
	for(i in 2:length(bulk_cl_hc_raw$cluster)){
		new_cluster[[i-1]] <- bulk_cl_hc_raw$cluster[[i]]
	}
	bulk_cl_hc <- list(reordering = match(avail_AD_genes[bulk_cl_hc_raw$reordering[-(1:3)]], genes),
                   cluster = new_cluster)
	# reload data based on filtered genes
	data_list <- load_data(genes)
	full_data_list <- list(data_list, genes)
}

saveRDS(list(full_data_list = full_data_list, bulk_cl_hc = bulk_cl_hc),
        sprintf('%s/%s_data.rds', output_dir, geneset))

