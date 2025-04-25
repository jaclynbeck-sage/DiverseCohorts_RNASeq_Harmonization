library(edgeR)
source("helper_functions.R")

metadata <- download_metadata()
counts <- download_rsem("syn64176431")

metadata <- subset(metadata, specimenID %in% colnames(counts))
counts <- counts[, metadata$specimenID]

stopifnot(all(duplicated(metadata$specimenID) == FALSE))

orig_size <- ncol(counts)

# Normalize data for PCA
data_orig <- DGEList(counts, samples = metadata, remove.zeros = TRUE)
data_orig <- normLibSizes(data_orig, method = "TMM")

counts_log <- log2(cpm(data_orig) + 0.5)

metadata <- validate_sex(metadata, counts_log)
metadata <- outlier_pca(metadata, counts_log)

metadata$lib_size <- colSums(counts)

ggplot(metadata, aes(x = tissue, y = lib_size, fill = tissue)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(size = 0.5) +
  theme_bw()

metadata <- validate_DV200(metadata)


# Save samples that passed QC

metadata$valid <- metadata$sex_valid & metadata$pca_valid & metadata$dv200_valid
print(table(metadata$tissue, metadata$valid))

metadata <- subset(metadata, valid == TRUE)
counts <- counts[, metadata$specimenID]

message(str_glue("{ncol(counts)} of {orig_size} samples passed QC."))

data_final <- DGEList(counts, samples = metadata, remove.zeros = TRUE)
data_final <- normLibSizes(data_final, method = "TMM")

saveRDS(data_final, file.path("data", "QC", "MSSM_qc.rds"))
