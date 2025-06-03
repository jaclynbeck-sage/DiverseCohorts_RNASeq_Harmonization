source("helper_functions.R")

raw_data <- readRDS(file.path("data", "QC", "Rush_qc.rds"))
cqn_data <- readRDS(file.path("data", "cqn", "Rush_cqn.rds"))

cqn_plt_data <- lapply(cqn_data, function(cqn_obj) {
  (cqn_obj$y + cqn_obj$offset) |>
  as.data.frame() |>
  tibble::rownames_to_column("ensembl_gene_id") |>
  tidyr::pivot_longer(cols = -ensembl_gene_id,
                      names_to = "specimenID",
                      values_to = "cqn_expr")
})

plts <- lapply(names(cqn_plt_data), function(tissue) {
  cqn_plt_tissue <- cqn_plt_data[[tissue]]

  ggplot(cqn_plt_tissue, aes(x = cqn_expr, color = specimenID)) +
    geom_density() +
    theme_bw() +
    theme(legend.position = "none") +
    ggtitle(tissue)
})

print(plts)


# Filter out highly correlated technical variables -----------------------------

numerical_vars <- raw_data$metadata |>
  select(specimenID, tissue, RIN, DV200, PMI) |>
  merge(raw_data$multiqc_stats) |>
  mutate(across(c(-specimenID, -tissue), as.numeric))

# Remove any variables where the values are the same for every sample
vars_remove <- numerical_vars |>
  summarize(across(where(is.numeric), ~ length(unique(.x)) == 1))

vars_remove <- c(colnames(vars_remove)[vars_remove == TRUE],
                 # Not really technical variables
                 "samtools_percent_mapped_X", "samtools_percent_mapped_Y",
                 # redundant with picard_PERCENT_DUPLICATION
                 "samtools_reads_duplicated_percent",
                 # Not enough variation
                 "samtools_average_quality",
                 # Not as useful for alignments to genome instead of transcriptome
                 "samtools_insert_size_average")

numerical_vars <- select(numerical_vars, -all_of(vars_remove))

numerical_vars <- lapply(unique(raw_data$metadata$tissue), function(tissue) {
  tech_mat <- numerical_vars[numerical_vars$tissue == tissue, ] |>
    select(-specimenID, -tissue) |>
    as.matrix()

  # Number of NA entries in each column
  na_vars <- numerical_vars[numerical_vars$tissue == tissue, ] |>
    summarize(across(everything(), ~sum(is.na(.x))))

  cor_mat <- cor(tech_mat, use = "na.or.complete")

  vars_remove <- remove_correlated_variables(cor_mat, na_vars)
  print(paste0(tissue, ": removing ", paste(vars_remove, collapse = ", ")))

  numerical_vars[numerical_vars$tissue == tissue, ] |>
    select(!any_of(vars_remove))
})

names(numerical_vars) <- unique(raw_data$metadata$tissue)

categorical_vars <- raw_data$metadata |>
  select(
    specimenID, tissue, rnaBatch, libraryBatch, sequencingBatch, isSampleExchange, cohort
    # sampleExchangeOrigin has 1:1 overlap with cohort so it is excluded
    # dataContributionGroup has 1:1 overlap with cohort and a near 1:1 overlap
    # with rnaBatch so it is excluded too.
  )

# Split by tissue, remove columns that are all NA, all have the same value, or
# have the same number of unique values as there are samples
categorical_vars <- lapply(unique(raw_data$metadata$tissue), function(tissue) {
  remove_vars <- categorical_vars[categorical_vars$tissue == tissue, ] |>
    summarize(
      across(
        -specimenID,
        ~ all(is.na(.x)) || length(unique(.x)) == 1 || length(unique(.x)) == length(.x)
      )
    )

  remove_vars <- names(remove_vars)[remove_vars == TRUE]

  categorical_vars[categorical_vars$tissue == tissue, ] |>
    select(!all_of(remove_vars)) |>
    mutate(across(-specimenID, factor))
})

names(categorical_vars) <- unique(raw_data$metadata$tissue)


for (tissue in unique(raw_data$metadata$tissue)) {
  # In case there is a random element anywhere in this process, set a seed for reproducibility
  set.seed(sum(utf8ToInt(paste("Rush", tissue, "regression"))))

  data_sub <- cqn_data[[tissue]]
  data_sub <- data_sub$y + data_sub$offset

  meta_sub <- merge(categorical_vars[[tissue]], numerical_vars[[tissue]])
  data_sub <- data_sub[, meta_sub$specimenID]

  stopifnot(all(colnames(data_sub) == meta_sub$specimenID))

  variables <- setdiff(colnames(meta_sub), c("specimenID", "tissue"))
  baseFormula <- "~1"

  mixed_vars <- intersect(variables,
                          c("cohort",
                            #"sequencingBatch", # For Rush there are only 2-3 levels of sequencing batch
                            "rnaBatch", "libraryBatch"))
  fixed_vars <- variables #setdiff(variables, mixed_vars)

  # All cohorts only belong to one dataContributionGroup each
  if (all(c("cohort", "dataContributionGroup") %in% variables)) {
    mixed_vars <- c(setdiff(mixed_vars, "cohort"),
                    "dataContributionGroup:cohort")
  }

  # For Rush (need to check other datasets), all RNA batches belong to one
  # dataContributionGroup each
  #if (all(c("rnaBatch", "dataContributionGroup") %in% variables)) {
  #  mixed_vars <- c(setdiff(mixed_vars, "rnaBatch"),
  #                  "dataContributionGroup:rnaBatch")
  #}

  mixed_vars <- paste0("(1|", mixed_vars, ")")

  meta_sub <- meta_sub |>
    # NA values for DV200 and RIN can be imputed as the median of their rna batch
    group_by(rnaBatch) |>
    mutate(across(any_of(c("DV200", "RIN")), ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x))) |>
    ungroup() |>
    mutate(
      # Other numeric variables with NAs (e.g. PMI) are set to the median of the whole column
      across(where(is.numeric), ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x)),
      across(where(is.numeric), scale),
      # Change NA values for factors to "missing or unknown"
      across(where(is.factor), ~ factor(ifelse(is.na(.x), "missing or unknown", as.character(.x))))
    )

  results <- mvIC::mvForwardStepwise(data_sub,
                                     baseFormula = baseFormula,
                                     data = meta_sub,
                                     variables = c(fixed_vars, mixed_vars))
  print(results$formula)

  saveRDS(results, file = file.path("data", "regression", paste0("Rush_", tissue, "_formulas.rds")))

  model_vars <- all.vars(results$formula)
  fixed_vars <- setdiff(model_vars,
                        c("cohort", "dataContributionGroup",
                          "libraryBatch", "rnaBatch"))#, "sequencingBatch"))
  mixed_vars <- setdiff(model_vars, fixed_vars)

  if (all(c("cohort", "dataContributionGroup") %in% model_vars)) {
    mixed_vars <- c(setdiff(mixed_vars, c("cohort", "dataContributionGroup")),
                    "dataContributionGroup/cohort")
  }

  if (all(c("rnaBatch", "dataContributionGroup") %in% model_vars)) {
    mixed_vars <- c(setdiff(mixed_vars, c("rnaBatch", "dataContributionGroup")),
                    "dataContributionGroup/rnaBatch")
  }

  mixed_vars <- paste0("(1 | ", mixed_vars, ")")

  formula_final <- paste("row ~", paste(c(fixed_vars, mixed_vars), collapse = " + "))

  print(paste0(tissue, ": ", formula_final))

  # LME4 can only do one gene at a time.
  # Hack to get row names back after mclapply, since I'm not sure the rows
  # are always returned in the same order
  genes <- rownames(data_sub)
  names(genes) <- genes

  n_cores <- parallel::detectCores()/2

  fits_lme <- lapply(genes, function(gene) { #parallel::mclapply(genes, function(gene) {
    row <- data_sub[gene, ]
    lme4::lmer(as.formula(formula_final), data = meta_sub)
  })#, mc.cores = n_cores)

  resid_lme <- t(sapply(fits_lme, residuals))
  resid_lme <- resid_lme[rownames(data_sub), ] # Ensure rows are in the original order
  colnames(resid_lme) <- colnames(data_sub)

  write.csv(resid_lme, file.path("data", "regression", str_glue("Rush_{tissue}_residuals.csv")),
            quote = FALSE)

  # Write original (un-scaled) metadata values
  meta_sub_orig <- merge(categorical_vars[[tissue]], numerical_vars[[tissue]])
  write.csv(meta_sub_orig[, c("specimenID", "tissue", model_vars)],
            file.path("data", "regression", str_glue("Rush_{tissue}_covariates.csv")),
            row.names = FALSE, quote = FALSE)

  rm(fits_lme)
  gc()
}

# TODO/notes:
# * consider libraryBatch or rnaBatch instead of sequencingBatch? a lot of the
#     sample swaps have NA for these values
# * caudate nucleus only has 2 sequencing batches -- doesn't need to be a mixed effect?
# * DLPFC -- 3 sequencing batches, same question
# * STG -- 2 sequencing batches, same question
# * Each dataContributionGroup is a linear combination of one or more cohorts.
#     Not sure if both variables are needed or just cohort. Could also be
#     specified as (1 | dataContributionGroup / cohort)
# * Need to set sampleExchangeOrigin to "not applicable" when it equals the
#     dataContributionSite (or MSSM/Columbia for NYGC)
# * Cohort and sampleExchangeOrigin are 1:1 so I don't think we need sampleExchangeOrigin
# * The samples with NA DV200 or RIN are mostly sample exchanges
# * Mayo and Emory are missing a large portion of PMI values (~ 42%)




