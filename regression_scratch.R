source("helper_functions.R")

raw_data <- readRDS(file.path("data", "QC", "Mayo_qc.rds"))
cqn_data <- readRDS(file.path("data", "cqn", "Mayo_cqn.rds"))

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

tech_vars <- raw_data$metadata |>
  select(specimenID, tissue, RIN, DV200, PMI) |>
  merge(raw_data$multiqc_stats) |>
  mutate(across(c(-specimenID, -tissue), as.numeric))

# Remove any variables where the values are the same for every sample
vars_remove <- tech_vars |>
  summarize(across(where(is.numeric), ~ length(unique(.x)) == 1))

vars_remove <- c(colnames(vars_remove)[vars_remove == TRUE],
                 # Not really technical variables
                 "samtools_percent_mapped_X", "samtools_percent_mapped_Y",
                 # redundant with picard_PERCENT_DUPLICATION
                 "samtools_reads_duplicated_percent")

tech_vars <- select(tech_vars, -all_of(vars_remove))

tech_vars <- lapply(unique(raw_data$metadata$tissue), function(tissue) {
  tech_mat <- tech_vars[tech_vars$tissue == tissue, ] |>
    select(-specimenID, -tissue) |>
    as.matrix()

  # Number of NA entries in each column
  na_vars <- tech_vars[tech_vars$tissue == tissue, ] |>
    summarize(across(everything(), ~sum(is.na(.x))))

  cor_mat <- cor(tech_mat, use = "na.or.complete")

  vars_remove <- remove_correlated_variables(cor_mat, na_vars)
  print(paste("Removing", paste(vars_remove, collapse = ", ")))

  tech_vars[tech_vars$tissue == tissue, ] |>
    select(!any_of(vars_remove))
})

names(tech_vars) <- unique(raw_data$metadata$tissue)


# Remove the warn and valid columns added by QC, non-binned ages, and a few
# other columns added during QC. Also remove columns that are technical
# covariates

bio_vars <- raw_data$metadata |>
  dplyr::rename(ageDeath_binned = ageBin) |>
  select(-contains("warn"), -contains("valid"), -ageDeathNumeric, -pmiNumeric,
         -pmiBin, -RINbin, -ageDeath, -pca_group, -individualID_AMPAD_1.0,
         -DV200, -PMI, -RIN, -rnaBatch, -libraryBatch)

# Split by tissue, remove columns that are all NA or all have the same value
bio_vars <- lapply(unique(raw_data$metadata$tissue), function(tissue) {
  remove_vars <- bio_vars[bio_vars$tissue == tissue, ] |>
    summarize(across(everything(), ~ all(is.na(.x)) || length(unique(.x)) == 1))

  remove_vars <- names(remove_vars)[remove_vars == TRUE]

  bio_vars[bio_vars$tissue == tissue, ] |>
    select(!all_of(remove_vars)) |>
    mutate(across(everything(), factor))
})

names(bio_vars) <- unique(raw_data$metadata$tissue)


tissue = "superior temporal gyrus"
data_sub = cqn_data[[tissue]]
data_sub = data_sub$y + data_sub$offset

meta_sub = merge(bio_vars[[tissue]], tech_vars[[tissue]])
data_sub = data_sub[, meta_sub$specimenID]

variables <- setdiff(colnames(meta_sub), c("individualID", "specimenID", "tissue"))
baseFormula <- "~ ADoutcome"

mixed_vars <- intersect(variables,
                        c("cohort", "sequencingBatch",
                          "sampleExchangeOrigin", "dataContributionGroup"))
fixed_vars <- setdiff(variables, mixed_vars)

mixed_vars <- paste0("(1|", mixed_vars, ")")

meta_sub <- meta_sub |>
  mutate(across(where(is.numeric), scale))

results <- mvIC::mvForwardStepwise(data_sub,
                                   baseFormula = baseFormula,
                                   data = meta_sub,
                                   variables = c(fixed_vars, mixed_vars))
