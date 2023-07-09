library(Seurat)
library(MIND)

rosmap_sc_data_dir <- '/gpfs/gibbs/pi/zhao/cs2629/ROSMAP/GeneExpression/snRNAseqPFC_BA10'
output_dir <- 'output'

# for deriving the priors in bMIND 
ROSMAP_sc <- readRDS(sprintf('%s/seurat_obj.rds', rosmap_sc_data_dir))
sc_counts <- GetAssayData(object = ROSMAP_sc, slot = "counts")

# obtain prior estimates using single cell data based on bMIND

# remove In cells because it raises errors in prior fitting
sc_In_ind <- which(ROSMAP_sc$celltype == 'In')
length(sc_In_ind)

# make reference data

# remove In cells from the data
# because it raises errors in prior fitting
ref <- as.matrix(sc_counts[,-sc_In_ind])

# make meta ref
ref_meta <- data.frame(sample = ROSMAP_sc$projid[-sc_In_ind],
                       cell_type = ROSMAP_sc$celltype[-sc_In_ind])

# compute prior
prior = get_prior(sc = ref, meta_sc = ref_meta)
# only 10911 genes have valid priors inferred from the single cell data
print(dim(prior$profile))
print(dim(prior$covariance))

match_prior_frac <- match(all_cts, colnames(prior$profile))
prior$profile <- prior$profile[,match_prior_frac]
prior$covariance <- prior$covariance[,match_prior_frac,match_prior_frac]


saveRDS('output/bMIND_prior.rds', prior)
