library(dplyr)
library(Matrix)

# -
# single cell data
# -
# Download single cell data from Synapse: https://www.synapse.org/#!Synapse:syn21261143

# extract count data
sc_counts <- readMM("filtered_count_matrix.mtx")
rownames(sc_counts) <- readLines("filtered_gene_row_names.txt")
mdata <- read.delim("filtered_column_metadata.txt")
celltype <- mdata$broad.cell.type
colnames(sc_counts) <- mdata$TAG

# subset random cells 
# CIBERSORTx does not recommend using a large number of cells
set.seed(1)
sub_size <- 5000
cell_ind <- sample.int(ncol(sc_counts), sub_size, replace = FALSE)

sub_counts <- as.matrix(sc_counts[, cell_ind])
colnames(sub_counts) <- celltype[cell_ind]

# cell type proportions before and after sampling
celltype %>% table %>% print
print(table(colnames(sub_counts)))

# manually append some End cell
End_cells <- which(celltype == 'End')
add_End_counts <- as.matrix(sc_counts[,End_cells[!End_cells %in% cell_ind]])
colnames(add_End_counts) <- rep('End', sum(!End_cells %in% cell_ind))
sub_counts <- cbind(sub_counts, add_End_counts)

# check if the max count is over 50
# otherwise data would be considered as log-transformed and got exponentiated by CIBERSORTx
print(max(sub_counts))

max_per_cell <- apply(sub_counts, 2, max)
# all the high expressions are from neurons
table(colnames(sub_counts), max_per_cell > 50) %>% print


sub_counts_df <- cbind(GeneSymbol = rownames(sc_counts),
                       sub_counts)


print(sprintf('Save single cell reference data with %i genes and %i cells',
              nrow(sub_counts_df), ncol(sub_counts_df) - 1))

write.table(sub_counts_df, sprintf('input/ROSMAP_sc_ref_all_group_%i.txt', 
                                   ncol(sub_counts)),
            quote=F, sep = '\t', row.names=F)


