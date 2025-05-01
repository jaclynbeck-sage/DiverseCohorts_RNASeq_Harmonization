library(edgeR)
library(cqn)
source("helper_functions.R")

data <- readRDS(file.path("data", "QC", "Mayo_qc.rds"))

# For the purposes of CQN, any diagnosis of "missing or unknown" is categorized
# as "Other" to groups of 1
data$samples$ADoutcome[data$samples$ADoutcome == "missing or unknown"] <- "Other"
data$samples$group <- factor(paste(data$samples$tissue, data$samples$ADoutcome))

# Filter to genes that are expressed in at least one group
genes_keep <- filterByExpr(data, group = data$samples$group)
data <- data[genes_keep, ]

gc_info <- read.csv(file.path("data", "gene_lengths_gc.csv")) |>
  mutate(ensembl_gene_id = str_replace(ensembl_gene_id, "\\.[0-9]+", "")) |>
  subset(ensembl_gene_id %in% rownames(data))

rownames(gc_info) <- gc_info$ensembl_gene_id
gc_info <- gc_info[rownames(data), ]

#data_split <- lapply(unique(data$samples$group), function(group) {
#  data[, data$samples$group == group]
#})

data_tissue <- lapply(unique(data$samples$tissue), function(tissue) {
  data[, data$samples$tissue == tissue]
})

cqn_data <- lapply(data_tissue, function(dt) {
  dt_cqn <- cqn(dt$counts,
                x = gc_info$percent_gc_content,
                lengths = gc_info$gene_length,
                lengthMethod = "smooth")
  dt_cqn$y + dt_cqn$offset
})

cqn_data <- do.call(cbind, cqn_data)

# Reverse the normalization operation done in cqn, which is:
# log2(counts + 1) - log2(lib_size / 10^6)
# Reversing, this would be:
# 2^(cqn_data + log2(lib_size / 10^6)) - 1
corr_counts <- sweep(cqn_data, 2,
                     log2(data$samples[colnames(cqn_data), "lib.size"] / 10^6),
                     "+")
corr_counts <- 2^corr_counts - 1
corr_counts[corr_counts < 0] <- 0
