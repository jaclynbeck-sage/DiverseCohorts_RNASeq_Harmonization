library(edgeR)
source("helper_functions.R")

metadata <- download_metadata()
counts <- download_rsem("syn64176419")

metadata <- subset(metadata, specimenID %in% colnames(counts))
counts <- counts[, metadata$specimenID]

stopifnot(all(duplicated(metadata$specimenID) == FALSE))

counts_log <- lognorm(counts)

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

data_final <- DGEList(counts, samples = metadata, remove.zeros = TRUE)
data_final <- normLibSizes(data_final, method = "TMM")

saveRDS(data_final, file.path("data", "QC", "Mayo_qc.rds"))
