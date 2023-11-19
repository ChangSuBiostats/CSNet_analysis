library(dplyr)
library(biomaRt) 

################
# Download ROSMAP data from Synapse
# bulk RNA-seq: https://www.synapse.org/Portal.html#!Synapse:syn3388564
# metadata: https://www.synapse.org/#!Synapse:syn3157322
###############

# directory that saves relevant files
ROSMAP_dir <- ''

# raw gene expression data in FPKM
ROSMAP_genes = read.table(file = paste(ROSMAP_dir, 'GeneExpression/bulk_brain_RNAseq/ROSMAP_RNAseq_FPKM_gene.tsv', sep = '/'), 
                          sep = '\t', header = TRUE)
ROSMAP_genes_df = as.data.frame(ROSMAP_genes)
# metadata
clinical_meta = read.csv(paste(ROSMAP_dir, "Metadata/ROSMAP_clinical.csv", sep = '/'))
bios_meta = read.csv(paste(ROSMAP_dir, "Metadata/ROSMAP_biospecimen_metadata.csv", sep = '/'))
assay_meta = read.csv(paste(ROSMAP_dir, "Metadata/ROSMAP_assay_RNAseq_metadata.csv", sep = '/'))


################
# Annotate gene names
# as only Ensemble IDs were provided in the original dataset
################

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = 104)
ens = ROSMAP_genes_df$gene_id
ens_lookup = gsub('\\.[0-9]*$', '', ens) # remove trailings

length(unique(ens))
length(unique(ens_lookup))

# set useCache=FALSE to deal with dplyr errors. See https://support.bioconductor.org/p/p132704/
# I tried updating Bioconductor. doesn't help
annotLookup = biomaRt::getBM(
    mart=ensembl,
    attributes=c("ensembl_gene_id", "gene_biotype", "external_gene_name", "start_position","end_position"),
    filter="ensembl_gene_id",
    values=ens_lookup,
    uniqueRows=TRUE,
    useCache = FALSE) 
# match it to original gene names,
# even if that introduces NA
annotLookup_matched <- annotLookup[match(ens_lookup, annotLookup$ensembl_gene_id),]
dim(annotLookup_matched)
nrow(annotLookup_matched) - sum(is.na(annotLookup_matched$external_gene_name))

saveRDS(annotLookup_matched, 'output/ROSMAP_annotLookup_new.rds')


################
# Clean ROSMAP bulk data
################

# check sample names, identify redundant/repeated samples
sample_names = gsub('X','', colnames(ROSMAP_genes_df))[-c(1,2)]
sample_names_no_trailing = gsub('.{2}$', '', sample_names)
length(sample_names_no_trailing)
## catch batch IDs in the sample names
sample_trailing = sapply(sample_names, function(x) strsplit(x, split = '_')[[1]][3])
table(sample_trailing, assay_meta$libraryBatch[match(sample_names_no_trailing, 
                                                     assay_meta$specimenID)])
## detect multiple samples for the same subject
table(sample_names_no_trailing)[table(sample_names_no_trailing) > 1]
sample_names[sample_names_no_trailing == '492_120515']
redundant_ID = which(sample_names_no_trailing == '492_120515')[c(2,3)] # this subject has 3 samples
## any redo?
sample_names[grepl('redo', sample_names)] # only 1. should be of reasonable quality (redone)

                         
# obtain covariates (including batch) for the bulk RNA-seq samples                        
## merge clinical_meta with assay_meta with specimen IDs
## match individual ID with specimen ID
individualID <- bios_meta$individualID[match(sample_names_no_trailing, bios_meta$specimenID)]
## # There are two specimens with empty individual ID
bios_meta$individualID[match(c('764_130520', '800_130701'), bios_meta$specimenID)] %>% as.character %>% print
sample_names_no_trailing[individualID == '']
sum(individualID == '')

## extract demographics from clinical_meta with individual ID
int_meta <- clinical_meta[match(individualID, clinical_meta$individualID),]
## get matched projid
int_meta[match(c('764_130520', '800_130701'), int_meta$specimenID), ] %>% print

## extract technical confounders from assay_meta with specimen ID
int_meta$RIN <- assay_meta$RIN[match(sample_names_no_trailing, assay_meta$specimenID)]
int_meta$libraryBatch <- assay_meta$libraryBatch[match(sample_names_no_trailing,
                                                       assay_meta$specimenID)]
                         
# subset 
ROSMAP_fpkm <- as.matrix(ROSMAP_genes_df[,-c(1,2)])
summary(ROSMAP_fpkm[1,])
# remove redudant samples
ROSMAP_fpkm_rm_redun <- ROSMAP_fpkm[, -redundant_ID]

# relabel specimen ID by projid
colnames(ROSMAP_fpkm_rm_redun) <- int_meta$projid
cov_meta_rm_redun <- cov_meta
                         
# split data by batch 
# batch 7 and 8 have strong batch effects compared to 1:6 as they were sequenced separately)
# Reference: 'These late samples were sequenced in batch 2 on plates 7 and 8.' from https://www.synapse.org/Portal.html#!Synapse:syn3388564
# and 'We used 82 samples from the ROSMAP project that were not included in the primary analysis.' from https://www.nature.com/articles/s41593-018-0154-9
ROSMAP_fpkm_rm_redun_train <- ROSMAP_fpkm_rm_redun[, !int_meta$libraryBatch %in% c(7,8)]
ROSMAP_fpkm_rm_redun_valid <- ROSMAP_fpkm_rm_redun[, int_meta$libraryBatch %in% c(7,8)]

head(ROSMAP_fpkm_rm_redun_train)

dim(ROSMAP_fpkm_rm_redun_train)
dim(ROSMAP_fpkm_rm_redun_valid)
ROSMAP_genes <- ROSMAP_genes_df[[1]]

saveRDS(list(train = ROSMAP_fpkm_rm_redun_train,
            valid = ROSMAP_fpkm_rm_redun_valid,
            gene = ROSMAP_genes), 'output/fpkm_unadj_by_batch.rds')

# save a version to run CIBERSORTx
# run deconvolution for the training samples
exp_m <- as.matrix(ROSMAP_fpkm_rm_redun_train)
exp_df <- cbind(GeneSymbol = ROSMAP_genes$external_gene_name,
                exp_m)
write.table(exp_df, sprintf('output/bulk_%s_%s.txt', set),
            row.names=F, quote=F, sep = '\t')


################
# Download ROSMAP snRNA-seq data from Synapse
# https://www.synapse.org/#!Synapse:syn21261143
###############

# load downloaded snRNA-seq data
sc_counts <- readMM("filtered_count_matrix.mtx")
rownames(sc_counts) <- readLines("filtered_gene_row_names.txt")
mdata <- read.delim("filtered_column_metadata.txt")


###############
# Generate cell-type-specific pseudo-bulk data
# to be used in compute_estimates.R
# for estimating single cell based co-expression
###############

unique_projid <- unique(mdata$projid)
for(ct in c('Ex', 'Oli', 'Mic', 'Ast')){
	count_m <- sc_counts[, mdata$broad.cell.type == ct]
	projid_ct <- mdata$projid[mdata$broad.cell.type == ct]
	syn_bulk <- matrix(0, nrow = nrow(count_m), ncol = length(unique_projid))
    rownames(syn_bulk) <- rownames(count_m)
    colnames(syn_bulk) <- unique_projid
    for(i in 1:length(unique_projid)){
        projid_ct_ind <- projid_ct == unique_projid[i] 
        if(sum(projid_ct_ind) > 1){
            syn_bulk[, i] <- rowSums(count_m[, projid_ct_ind])
        }else if(sum(projid_ct_ind) == 1){
            syn_bulk[, i] <- count_m[, projid_ct_ind]
        }else if(sum(projid_ct_ind) == 0){
            syn_bulk[, i] <- rep(0, nrow(syn_bulk))
        }
    }
    seq_depths <- apply(syn_bulk, 2, sum)

    # rescale by normalizing depths
    # to have the same magnitude: median(seq_depths)
    syn_bulk <- scale(syn_bulk, scale = seq_depths, center = F) * median(seq_depths)
    saveRDS(sprintf('output/pseudo_bulk/bulk_scRNAseq_%s_nor_pseudo_bulk.rds', ct))
}
