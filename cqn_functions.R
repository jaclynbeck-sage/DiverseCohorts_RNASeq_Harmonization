# This script is set up so that sections of it will be embedded in the Quarto
# CQN notebooks. All this code is here instead of directly in the notebooks
# because there are 4 datasets / 4 notebooks, and moving the code here allows me
# to change something in one place and have it propagate to all 4 notebooks
# instead of having to change it 4 times.

# ---- include-libraries ----

library(edgeR)
library(cqn)
library(stringr)
library(sageRNAUtils)
library(ggplot2)
library(viridis)
library(patchwork)
library(synapser)

configs <- config::get(file = "config.yml")

# Defaults that should be edited in the qmd file
dataset <- NULL
upload_to_synapse <- FALSE


# ---- load-data ----

synLogin()

gene_file <- synGet(configs$download$gene_metadata_synid,
                    downloadLocation = "downloads")
gene_info <- read.csv(gene_file$path)

data <- readRDS(file.path("data", "QC", paste0(dataset, "_qc.rds")))
data <- DGEList(round(data$counts), samples = data$metadata)

# For the purposes of filterByExpr, any diagnosis of "missing or unknown" is
# categorized as "Other" to avoid groups of 1
data$samples$ADoutcome[data$samples$ADoutcome == "missing or unknown"] <- "Other"
data$samples$group <- data$samples$ADoutcome

# Filter to genes that are expressed in at least one group, for each tissue separately
data_tissue <- lapply(unique(data$samples$tissue), function(tissue) {
  data_sub <- data[, data$samples$tissue == tissue]
  genes_keep <- edgeR::filterByExpr(data_sub, group = data_sub$samples$group)
  data_sub[genes_keep, ]
})

names(data_tissue) <- unique(data$samples$tissue)
rownames(gene_info) <- gene_info$ensembl_gene_id


# ---- run-cqn ----

cqn_data <- lapply(data_tissue, function(dt) {
  # There seems to be a random element to cqn() so we set a seed for reproducibility
  set.seed(sageRNAUtils::string_to_seed(
    paste(dataset, unique(dt$samples$tissue), "cqn")
  ))

  gene_info_sub <- gene_info[rownames(dt$counts), ]

  stopifnot(all(rownames(dt$counts) == gene_info_sub$ensembl_gene_id))

  output <- cqn(dt$counts,
                x = gene_info_sub$percent_gc_content,
                lengths = gene_info_sub$gene_length,
                lengthMethod = "smooth")
  output[c("y", "offset", "glm.offset")]
})

saveRDS(cqn_data, file.path("data", "cqn", paste0(dataset, "_cqn.rds")))

# Write each tissue's CQN data to a CSV file and optionally upload to Synapse
for (tissue in names(cqn_data)) {
  cqn_file <- file.path("data", "cqn", str_glue("{dataset}_{tissue}_cqn.csv"))

  write.csv(cqn_data[[tissue]]$y + cqn_data[[tissue]]$offset,
            cqn_file,
            quote = FALSE)

  if (upload_to_synapse) {
    provenance <- c(
      synFindEntityId(name = paste0(dataset, "_counts_filtered.csv"),
                      parent = configs$upload$counts_folder_synid),
      configs$download$gene_metadata_synid
    )

    github <- c(
      paste0(
        "https://github.com/jaclynbeck-sage/DiverseCohorts_RNASeq_Harmonization/",
        "blob/main/03_", dataset, "_CQN.qmd"
      ),
      paste0(
        "https://github.com/jaclynbeck-sage/DiverseCohorts_RNASeq_Harmonization/",
        "blob/main/cqn_functions.R"
      )
    )

    dataset_folder <- Folder(dataset, parent = configs$upload$cqn_folder_synid)
    dataset_folder <- synStore(dataset_folder)

    syn_file <- File(cqn_file, parent = dataset_folder)

    syn_file <- synStore(
      syn_file,
      forceVersion = FALSE,
      used = provenance,
      executed = github
    )
  }
}


# ---- plot-cqn ----

cqn_plt_data <- lapply(cqn_data, function(cqn_obj) {
  (cqn_obj$y + cqn_obj$offset) |>
    as.data.frame() |>
    tibble::rownames_to_column("ensembl_gene_id") |>
    tidyr::pivot_longer(cols = -ensembl_gene_id,
                        names_to = "specimenID",
                        values_to = "cqn_expr")
})

raw_plt_data <- lapply(cqn_data, function(cqn_obj) {
  cqn_obj$y |> # Un-adjusted log2-cpm data
    as.data.frame() |>
    tibble::rownames_to_column("ensembl_gene_id") |>
    tidyr::pivot_longer(cols = -ensembl_gene_id,
                        names_to = "specimenID",
                        values_to = "raw_expr")
})

plts <- lapply(names(cqn_plt_data), function(tissue) {
  cqn_plt_tissue <- cqn_plt_data[[tissue]]
  raw_plt_tissue <- raw_plt_data[[tissue]]

  plt1 <- ggplot(raw_plt_tissue, aes(x = raw_expr, color = specimenID)) +
    geom_density() +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = "Raw log2-CPM data", x = "Raw log2 expression", y = "Density") +
    scale_color_viridis(discrete = TRUE, begin = 0, end = 0.95)

  plt2 <- ggplot(cqn_plt_tissue, aes(x = cqn_expr, color = specimenID)) +
    geom_density() +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = "CQN", x = "CQN expression", y = "Density") +
    scale_color_viridis(discrete = TRUE, begin = 0, end = 0.95)

  print(plt1 + plt2 + plot_annotation(title = str_to_title(tissue)))
})
