library(dplyr)
library(stringr)
library(matrixStats)
library(sageRNAUtils)
library(variancePartition)

source("helper_functions.R")

dataset <- "Rush"

raw_data <- readRDS(file.path("data", "QC", str_glue("{dataset}_qc.rds")))
cqn_data <- readRDS(file.path("data", "cqn", str_glue("{dataset}_cqn.rds")))
clean_covariates <- readRDS(file.path("data", "regression",
                                      str_glue("{dataset}_cleaned_covariates.rds")))

n_cores <- parallel::detectCores() / 2


# Look at sources of variation in each tissue ----------------------------------

#fqc_tmp <- fastqc_data$basic_statistics |>
#  group_by(specimenID, tissue) |>
#  summarize(percent_gc_content_fastqc = mean(percent_gc_content), .groups = "drop")

#meta_merged <- metadata |>
#  merge(multiqc_stats) |>
#  merge(fqc_tmp)

#rownames(meta_merged) <- meta_merged$specimenID

#form <- ~ (1 | sex) + (1 | race) + (1 | isHispanic) + (1 | ADoutcome) +
#  (1 | rnaBatch) + RIN + DV200 + picard_PERCENT_DUPLICATION +
#  percent_gc_content_fastqc + rsem_alignable_percent

plot_sources_of_variance(cqn_data, clean_covariates)

formulas <- list(
  "dorsolateral prefrontal cortex" = ~ DV200 + rsem_uniquely_aligned_percent +
    rsem_alignable_percent + picard_PERCENT_DUPLICATION + (1 | sequencingBatch) +
    (1 | cohort),
  "test" = ~ RIN + rsem_uniquely_aligned_percent + percent_gc_content_fastq +
    rsem_alignable_percent + (1 | libraryBatch) + (1 | cohort)
)


# Regress out technical variation from each tissue -----------------------------

for (tissue in names(clean_covariates)) {
  # In case there is a random element anywhere in this process, set a seed for reproducibility
  set.seed(sageRNAUtils::string_to_seed(paste(dataset, tissue, "regression")))

  data_sub <- cqn_data[[tissue]]
  data_sub <- data_sub$y + data_sub$offset

  meta_sub <- clean_covariates[[tissue]]
  data_sub <- data_sub[, meta_sub$specimenID]

  rownames(meta_sub) <- meta_sub$specimenID

  stopifnot(all(colnames(data_sub) == meta_sub$specimenID))

  mvIC_results <- readRDS(file.path("data", "regression",
                                    str_glue("{dataset}_{tissue}_mvIC_results.rds")))

  # Remove biological variables, which we aren't regressing out, and cohort,
  # which doesn't contribute much to variance in comparison to batch but may
  # have been added anyway
  model_vars <- all.vars(mvIC_results$formula) |>
    setdiff(c("ADoutcome", "sex", "race", "isHispanic", "ageDeath", "apoeGenotype", "cohort"))

  fixed_vars <- setdiff(model_vars,
                        c("libraryBatch", "rnaBatch", "sequencingBatch"))

  mixed_vars <- setdiff(model_vars, fixed_vars)

  if (length(mixed_vars) > 0) {
    mixed_vars <- paste0("(1 | ", mixed_vars, ")")
  }

  form <- paste("~", paste(c(fixed_vars, mixed_vars), collapse = " + "))
  print(paste0(tissue, ": ", form))

  fit <- fitVarPartModel(data_sub, form, meta_sub,
                         BPPARAM = MulticoreParam(n_cores))

  resid_tissue <- residuals(fit)
  rownames(resid_tissue) <- rownames(data_sub)

  plot_resid_qq(resid_tissue)

  write.csv(resid_tissue,
            file.path("data", "regression", str_glue("{dataset}_{tissue}_residuals.csv")),
            quote = FALSE)

  # Write metadata values
  model_vars <- all.vars(as.formula(form))
  write.csv(meta_sub[, c("specimenID", model_vars)],
            file.path("data", "regression", str_glue("{dataset}_{tissue}_covariates.csv")),
            row.names = FALSE, quote = FALSE)

  rm(fit)
  gc()
}

# TODO/notes:
# * caudate nucleus only has 2 sequencing batches -- doesn't need to be a mixed effect?
# * DLPFC -- 3 sequencing batches, same question
# * STG -- 2 sequencing batches, same question
# * Need to set sampleExchangeOrigin to "not applicable" when it equals the
#     dataContributionSite (or MSSM/Columbia for NYGC)
# * Cohort and sampleExchangeOrigin are 1:1 so I don't think we need sampleExchangeOrigin
# * The samples with NA DV200 or RIN are mostly sample exchanges
# * Mayo and Emory are missing a large portion of PMI values (~ 42%)




