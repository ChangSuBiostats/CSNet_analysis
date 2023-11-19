library(MIND)
library(Matrix)

rosmap_sc_data_dir <- '/gpfs/gibbs/pi/zhao/cs2629/ROSMAP/GeneExpression/snRNAseqPFC_BA10'
output_dir <- 'output'

# load single cell data
sc_counts <- readMM("filtered_count_matrix.mtx")
rownames(sc_counts) <- readLines("filtered_gene_row_names.txt")
mdata <- read.delim("filtered_column_metadata.txt")
celltype <- mdata$broad.cell.type
projid <- mdata$projid
colnames(sc_counts) <- mdata$TAG

# obtain prior estimates using single cell data based on bMIND

# remove In cells because it raises errors in prior fitting
sc_In_ind <- which(celltype == 'In')
length(sc_In_ind)

# make reference data

# remove In cells from the data
# because it raises errors in prior fitting
ref <- as.matrix(sc_counts[,-sc_In_ind])

# make meta ref
ref_meta <- data.frame(sample = projid[-sc_In_ind],
                       cell_type = celltype[-sc_In_ind])

# compute prior
prior = get_prior(sc = ref, meta_sc = ref_meta)
# only 10911 genes have valid priors inferred from the single cell data
print(dim(prior$profile))
print(dim(prior$covariance))

match_prior_frac <- match(all_cts, colnames(prior$profile))
prior$profile <- prior$profile[,match_prior_frac]
prior$covariance <- prior$covariance[,match_prior_frac,match_prior_frac]


saveRDS('output/bMIND_prior.rds', prior)
