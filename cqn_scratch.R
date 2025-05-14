library(edgeR)
library(cqn)
source("helper_functions.R")

data <- readRDS(file.path("data", "QC", "Rush_qc.rds"))

data <- DGEList(data$counts, samples = data$metadata)

# For the purposes of CQN, any diagnosis of "missing or unknown" is categorized
# as "Other" to avoid groups of 1
data$samples$ADoutcome[data$samples$ADoutcome == "missing or unknown"] <- "Other"
data$samples$group <- factor(paste(data$samples$tissue, data$samples$ADoutcome))


# Filter to genes that are expressed in at least one group
genes_keep <- filterByExpr(data, group = data$samples$group)
data <- data[genes_keep, ]

gc_info <- read.csv(file.path("data", "gene_metadata.csv"))

rownames(gc_info) <- gc_info$ensembl_gene_id
gc_info <- gc_info[rownames(data), ]

#data_split <- lapply(unique(data$samples$group), function(group) {
#  data[, data$samples$group == group]
#})

data_tissue <- lapply(unique(data$samples$tissue), function(tissue) {
  data[, data$samples$tissue == tissue]
})

cqn_data <- lapply(data_tissue, function(dt) {
  cqn(dt$counts,
      x = gc_info$percent_gc_content,
      lengths = gc_info$gene_length,
      lengthMethod = "smooth")
})

cqn_matrix <- do.call(cbind, lapply(cqn_data, function(cqn_obj) {
  cqn_obj$y + cqn_obj$offset
}))

cqn_matrix <- cqn_matrix[, data$samples$specimenID]

write.csv(cqn_matrix, file.path("data", "cqn", "Rush_cqn.csv"), quote = FALSE)

# Reverse the normalization operation done in cqn
corr_counts <- cqn_to_counts(cqn_matrix, data$samples$lib.size)
