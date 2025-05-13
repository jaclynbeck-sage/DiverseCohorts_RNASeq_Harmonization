library(edgeR)
source("helper_functions.R")

# TODO counts matrix is missing 26 samples, all sample exchanges
metadata <- download_metadata()
counts <- download_rsem("syn64176419")

fastqc_data <- readRDS(file.path("data", "QC", "Mayo_Emory_fastqc_stats.rds"))

multiqc_stats <- readRDS(file.path("data", "QC", "Mayo_Emory_multiqc_stats.rds"))

gene_info <- read.csv(file.path("data", "gene_lengths_gc.csv"))

metadata <- subset(metadata, specimenID %in% colnames(counts))
counts <- counts[, metadata$specimenID]

fastqc_data <- lapply(fastqc_data, function(df) {
  merge(dplyr::select(metadata, specimenID, tissue), df)
})

multiqc_stats <- merge(dplyr::select(metadata, specimenID, tissue), multiqc_stats)

stopifnot(length(unique(metadata$specimenID)) == nrow(metadata))

orig_size <- ncol(counts)

counts_log <- simple_lognorm(counts)

metadata <- validate_fastqc(metadata, fastqc_data)
metadata <- validate_sex(metadata, counts_log)
metadata <- outlier_pca(metadata, counts_log, gene_info)

metadata$lib_size <- colSums(counts)

ggplot(metadata, aes(x = tissue, y = lib_size, fill = tissue)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(size = 0.5) +
  theme_bw()

metadata <- validate_DV200(metadata)


# Save samples that passed QC

metadata$valid <- metadata$phred_score_valid & metadata$sex_valid &
  metadata$pca_valid & metadata$dv200_valid
# TODO warnings

print(table(metadata$tissue, metadata$valid))

metadata <- subset(metadata, valid == TRUE)
counts <- counts[, metadata$specimenID]

message(str_glue("{ncol(counts)} of {orig_size} samples passed QC."))

data_final <- DGEList(counts, samples = metadata, remove.zeros = TRUE)
#data_final <- normLibSizes(data_final, method = "TMM")

saveRDS(data_final, file.path("data", "QC", "Mayo_qc.rds"))
