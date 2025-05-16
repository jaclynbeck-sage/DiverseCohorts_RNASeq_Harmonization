library(edgeR)
library(cqn)
source("helper_functions.R")

gc_info <- read.csv(file.path("data", "gene_metadata.csv"))

data <- readRDS(file.path("data", "QC", "Mayo_qc.rds"))
data <- DGEList(round(data$counts), samples = data$metadata)

# For the purposes of CQN, any diagnosis of "missing or unknown" is categorized
# as "Other" to avoid groups of 1
data$samples$ADoutcome[data$samples$ADoutcome == "missing or unknown"] <- "Other"
data$samples$group <- data$samples$ADoutcome

# Filter to genes that are expressed in at least one group
data_tissue <- lapply(unique(data$samples$tissue), function(tissue) {
  data_sub <- data[, data$samples$tissue == tissue]
  genes_keep <- filterByExpr(data_sub, group = data_sub$samples$group)
  data_sub[genes_keep, ]
})

names(data_tissue) <- unique(data$samples$tissue)

rownames(gc_info) <- gc_info$ensembl_gene_id

cqn_data <- lapply(data_tissue, function(dt) {
  gc_info_sub <- gc_info[rownames(dt$counts), ]
  output <- cqn(dt$counts,
                x = gc_info_sub$percent_gc_content,
                lengths = gc_info_sub$gene_length,
                lengthMethod = "smooth")
  output[c("y", "offset", "glm.offset")]
})

saveRDS(cqn_data, file.path("data", "cqn", "Mayo_cqn.rds"))
