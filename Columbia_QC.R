library(edgeR)
source("helper_functions.R")

metadata <- download_metadata()
counts <- download_rsem("syn64289221")

fastqc_data <- readRDS(file.path("data", "QC", "Columbia_fastqc_stats.rds"))
fastqc_stats <- fastqc_data$fastq_summary
gc_distribution <- fastqc_data$gc_distribution

multiqc_stats <- readRDS(file.path("data", "QC", "Columbia_multiqc_stats.rds"))

gene_info <- read.csv(file.path("data", "gene_lengths_gc.csv"))

metadata <- subset(metadata, specimenID %in% colnames(counts))
counts <- counts[, metadata$specimenID]
fastqc_stats <- subset(fastqc_stats, specimenID %in% metadata$specimenID)
multiqc_stats <- subset(multiqc_stats, specimenID %in% metadata$specimenID)

stopifnot(all(duplicated(metadata$specimenID) == FALSE))

orig_size <- ncol(counts)

counts_log <- lognorm(counts)

metadata <- validate_sex(metadata, counts_log)
metadata <- outlier_pca(metadata, counts_log, gene_info)

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

saveRDS(data_final, file.path("data", "QC", "Columbia_qc.rds"))
