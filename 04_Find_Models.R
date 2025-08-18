# This process can take a long time, up to several hours depending on the number
# of cores the machine has, so it is not run as a notebook. Rather it is run
# separately as a script and the results are examined in a notebook in step 05.
library(dplyr)
library(stringr)
library(sageRNAUtils)
library(mvIC)
library(variancePartition)

source("helper_functions.R")

dataset <- "Mayo_Emory"
n_cores <- parallel::detectCores() - 2

raw_data <- readRDS(file.path("data", "QC", str_glue("{dataset}_qc.rds")))
cqn_data <- readRDS(file.path("data", "cqn", str_glue("{dataset}_cqn.rds")))

bio_vars <- c("ADoutcome", "ageDeath", "apoeGenotype", "isHispanic", "race",
              "sex")


# Clean up and filter technical variables --------------------------------------

numerical_vars <- raw_data$metadata |>
  dplyr::select(specimenID, tissue, RIN, DV200, PMI) |>

  # Ensure all data is numeric
  mutate(across(c(-specimenID, -tissue), as.numeric))

categorical_vars <- raw_data$metadata |>
  dplyr::select(
    specimenID, tissue, all_of(bio_vars),
    # Technical variables
    rnaBatch, libraryBatch, sequencingBatch, cohort
  ) |>
  mutate(ageDeath = bin_ages(ageDeath))

all_covariates <- merge(numerical_vars, categorical_vars)

if (dataset == "MSSM") {
  all_covariates <- all_covariates |>
    merge(select(raw_data$metadata, specimenID, dataContributionGroup))
}

clean_covariates <- lapply(unique(raw_data$metadata$tissue), function(tissue) {
  # Impute NA values and scale numeric covariates.
  covar <- all_covariates[all_covariates$tissue == tissue, ] |>

    # DV200 and RIN can be imputed as the median of their rna batch. This is
    # done before removing unusable covariates in case RNA batch gets removed
    # (which happens with Columbia data).
    group_by(rnaBatch) |>
    mutate(across(any_of(c("DV200", "RIN")),
                  ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x))) |>
    ungroup() |>

    # Remove covariates that have all NAs or only 1 unique value in the column
    remove_unusable_covariates(always_keep = c("specimenID", "tissue")) |>

    mutate(
      # Other numeric variables with NAs (e.g. PMI) are set to the median of the whole column
      across(c(where(is.numeric), -any_of(c("DV200", "RIN"))),
             ~ ifelse(is.na(.x), median(.x, na.rm = TRUE), .x)),

      # Scale
      across(where(is.numeric), ~ as.numeric(scale(.x))),

      # Change NA values for factors to "missing or unknown"
      across(c(where(is.character), where(is.factor), -specimenID),
             ~ factor(ifelse(is.na(.x), "missing or unknown", as.character(.x))))
    ) |>

    # Convert from tibble -> data.frame
    as.data.frame()

  # Mayo/Emory, Columbia, and MSSM:
  # "race" and "isHispanic" have a near-1:1 relationship, so we combine them
  # into one variable (race).
  #if (dataset %in% c("Mayo_Emory", "Columbia", "MSSM")) {
    covar <- covar |>
      mutate(
        isHispanic = case_match(
          as.character(isHispanic),
          "True" ~ "Hispanic",
          "False" ~ "Not Hispanic",
          .default = "missing or unknown"
        ),
        race = paste(race, "/", isHispanic)
      ) |>
      select(-isHispanic)
  #}

  # Mayo/Emory only:
  # - rnaBatch, libraryBatch, and sequencingBatch are all the same values for
  #   all tissues so we remove two of the three variables.
  # - race is split very unevenly across cohorts, with some cohorts having almost
  #   all samples from one race, so cohort is removed as a covariate to avoid
  #   interaction with race
  # - Remove PMI because of too many missing values
  if (dataset == "Mayo_Emory") {
     covar <- covar |>
       select(-rnaBatch, -libraryBatch, -cohort, -PMI)
  }

  # Columbia only:
  # - remove libraryBatch as a variable: race is spread very
  #   unevenly across batches, with some batches having samples from only one
  #   race.
  # - Some of Columbia's batches have a 1:1 or 1:many relationship:
  #   CN: 1 rnaBatch : 1 sequencingBatch, only 2 rnaBatches / sequencingBatches
  #       total so both variables are removed.
  #   DLPFC: rnaBatch is generally 1:1 with libraryBatch except for libraryBatch
  #          "NYGC_Columbia_B02", which is spread across 5 rnaBatches.
  #          sequencingBatch is generally 1:1 with libraryBatch except for
  #          libraryBatch "NYGC_Columbia_B02" (same as above), which is spread
  #          across 3 sequencingBatches. The near 1:1 relationship means that
  #          race is spread unevenly across both these variables, so they are
  #          removed.
  #   TP: Only 1 rnaBatch and sequencingBatch, the columns were already removed
  #       above by remove_unusable_covariates().
  # - cohort is removed as a variable: there are only 2 main cohorts, 3 smaller
  #   cohorts with 1-3 samples per tissue. The small cohorts each only have
  #   samples from one race.
  if (dataset == "Columbia") {
    covar <- covar |>
      select(-cohort, -libraryBatch,
             -any_of(c("rnaBatch", "sequencingBatch"))) # These were removed already for the TP
  }

  # MSSM only:
  # - Remove cohort as a variable since there is only the main MSSM cohort plus
  #   5-6 sample swaps per tissue. Cohort is replaced with dataContributionGroup,
  #   which is added to the data frame outside this loop, for MSSM only.
  # - Remove rnaBatch as a variable: all of the race = White samples from MSSM
  #   are in one rnaBatch ("NYGC_MSSM_B03"), which also happens to be the batch
  #   with very low RIN + very high DV200. This batch has only 3 other non-White
  #   samples. The only other White samples in other batches are from Rush. This
  #   level of co-confounding makes it extremely difficult to figure out what to
  #   regress out.
  # - Remove sequencingBatch as a variable for the same reason. libraryBatches
  #   are, thankfully, balanced, but for some reason rna and sequencingBatches
  #   were divided by race.
  if (dataset == "MSSM") {
    covar <- covar |>
      select(-cohort, -rnaBatch)
  }

  # Rush only:
  # - There are mostly 2 libraryBatches to 1 rnaBatch, with very few
  #   libraryBatches that have more than one rnaBatch. sequencingBatch is
  #   almost exactly 1:many with rnaBatch. I've chosen to use libraryBatch for
  #   the CN and STG as it is less correlated with race. There is an extremely uneven distribution of race across
  #   some rna and libraryBatches in the DLPFC, so I've decided to use
  #   sequencingBatch for that tissue instead even though there are only 3 batches.
  # - Remove cohort as a variable: the MARS cohort has only samples that are
  #   Black / African American, and some cohorts only have 1-2 samples total.
  if (dataset == "Rush") {
    covar <- covar |>
      select(-cohort, -rnaBatch)

    if (tissue == "dorsolateral prefrontal cortex") {
      covar <- covar |>
        select(-libraryBatch)
    } else {
      covar <- covar |>
        select(-sequencingBatch)
    }
  }

  return(covar)
})

names(clean_covariates) <- unique(raw_data$metadata$tissue)

saveRDS(clean_covariates,
        file.path("data", "regression",
                  str_glue("{dataset}_cleaned_covariates.rds")))


# Calculate % of variance explained by each covariate --------------------------

# TODO for Rush DLPFC, sequencing batch is a linear combination of rnaBatch or
# libraryBatch or rnaBatch + libraryBatch so it doesn't add extra information.
# libraryBatch is *generally* a linear combination of rnaBatch, but rnaBatch
# is more correlated with other variables than libraryBatch, unclear which to use
var_part_list <- lapply(clean_covariates, function(covariates) {
  tissue <- unique(covariates$tissue) |> as.character()

  # In case there is a random element anywhere in this process, set a seed for reproducibility
  set.seed(sageRNAUtils::string_to_seed(paste(dataset, tissue, "varpart modeling")))

  data_tissue <- cqn_data[[tissue]]$y + cqn_data[[tissue]]$offset
  data_tissue <- data_tissue[, covariates$specimenID]

  # Do this on the top 5000 variable genes only so it doesn't take too long to run
  var_genes <- rowVars(data_tissue) |> sort(decreasing = TRUE) |> names()

  data_tissue <- data_tissue[var_genes[1:5000], ]

  # Remove specimenID and tissue from the possible covariates
  covariates <- covariates |>
    tibble::column_to_rownames("specimenID") |>
    select(-tissue)

  # Just to see a printout of what this function would remove. Doesn't actually
  # remove anything.
  tmp <- covariates |> sageRNAUtils::remove_correlated_covariates()

  covariates_tech <- covariates |>
    select(-any_of(bio_vars))

  # All categorical values have to be mixed variables for fitExtractVarPartModel
  mixed_vars <- covariates |>
    select(where(is.character), where(is.factor)) |>
    colnames()

  mixed_vars_tech <- setdiff(mixed_vars, bio_vars)

  fixed_vars <- setdiff(colnames(covariates), mixed_vars)

  if (length(mixed_vars) > 0) {
    mixed_vars <- paste0("(1 | ", mixed_vars, ")")
  }

  if (length(mixed_vars_tech) > 0) {
    mixed_vars_tech <- paste0("(1 | ", mixed_vars_tech, ")")
  }

  form_all <- paste("~", paste(c(fixed_vars, mixed_vars), collapse = " + "))
  form_tech <- paste("~", paste(c(fixed_vars, mixed_vars_tech), collapse = " + "))

  var_part1 <- fitExtractVarPartModel(data_tissue, form_all, covariates,
                                      BPPARAM = MulticoreParam(n_cores)) |>
    sortCols()

  var_part2 <- fitExtractVarPartModel(data_tissue, form_tech, covariates,
                                      BPPARAM = MulticoreParam(n_cores)) |>
    sortCols()

  plt <- plotVarPart(var_part1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(title = paste("Sources of variance:", tissue),
         subtitle = "All variables, sorted by median")
  print(plt)

  plt <- plotVarPart(sortCols(var_part1, FUN = mean)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(title = paste("Sources of variance:", tissue),
         subtitle = "All variables, sorted by mean")
  print(plt)

  plt <- plotVarPart(var_part2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(title = paste("Sources of variance:", tissue),
         subtitle = "Technical variables, sorted by median")
  print(plt)

  plt <- plotVarPart(sortCols(var_part2, FUN = mean)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(title = paste("Sources of variance:", tissue),
         subtitle = "Technical variables, sorted by mean")
  print(plt)

  return(list(var_part_all = var_part1,
              var_part_tech = var_part2))
})

vp_summary <- function(var_part) {
  var_part |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    select(-Residuals) |>
    tidyr::pivot_longer(-gene, names_to = "covariate", values_to = "value") |>
    group_by(covariate) |>
    summarize(mean_pct = mean(value),
              median_pct = median(value)) |>
    mutate(rank_mean = rank(-mean_pct),
           rank_median = rank(-median_pct),
           total_rank = (rank_mean + rank_median) / 2) |>
    arrange(desc(mean_pct))
}


for (tissue in names(clean_covariates)) {
  vp <- vp_summary(var_part_list[[tissue]]$var_part_all)
  print(vp)

  covariates <- clean_covariates[[tissue]] |>
    tibble::column_to_rownames("specimenID") |>
    select(-tissue)

  form_cc <- paste("~", paste0(colnames(covariates), collapse = " + "))
  cc <- canCorPairs(form_cc, covariates, showWarnings = FALSE)
  cc[upper.tri(cc, diag = TRUE)] <- NA

  print(cc)

  cor_pairs <- cc |>
    as.data.frame() |>
    tibble::rownames_to_column("var1") |>
    tidyr::pivot_longer(-var1, names_to = "var2", values_to = "cor_value") |>
    subset(cor_value > 0.7) |>
    arrange(desc(cor_value)) |>
    rowwise() |>
    mutate(
      var1_rank = vp$rank_mean[vp$covariate == var1],
      var2_rank = vp$rank_mean[vp$covariate == var2],
      lowest_pct = which.max(c(var1_rank, var2_rank)),
      var_remove = c(var1, var2)[lowest_pct]
    ) |>
    ungroup()

  print(cor_pairs)

  to_remove <- c()
  for (N in 1:nrow(cor_pairs)) {
    pair <- c(cor_pairs$var1[N], cor_pairs$var2[N])

    # Only remove one of the two variables if neither has been removed yet
    if (!any(pair %in% to_remove)) {
      to_remove <- c(to_remove, cor_pairs$var_remove[N])
    }
  }

  print(paste("Removing", length(to_remove), "variables:",
              paste(to_remove, collapse = ", ")))

  vars_keep <- vp |>
    arrange(desc(mean_pct)) |>
    subset(mean_pct > 0.01 & !(covariate %in% c(to_remove, bio_vars))) |>
    pull(covariate)

  print(paste(tissue, "- recommended variables:",
              paste(vars_keep, collapse = ", ")))

  saveRDS(list(recommended_vars = vars_keep,
               var_part = var_part_list[[tissue]],
               vp_summary_all = vp,
               vp_summary_tech = vp_summary(var_part_list[[tissue]]$var_part_tech),
               correlations = cc,
               cor_pairs = cor_pairs),
          file.path("data", "regression",
                    str_glue("{dataset}_{tissue}_var_part_results.rds")))
}


# Find a potential model for each tissue using mvIC ----------------------------

for (tissue in unique(raw_data$metadata$tissue)) {
  # In case there is a random element anywhere in this process, set a seed for reproducibility
  set.seed(sageRNAUtils::string_to_seed(paste(dataset, tissue, "mvIC modeling")))

  data_sub <- cqn_data[[tissue]]
  data_sub <- data_sub$y + data_sub$offset

  meta_sub <- clean_covariates[[tissue]]
  data_sub <- data_sub[, meta_sub$specimenID]

  stopifnot(all(colnames(data_sub) == meta_sub$specimenID))

  variables <- setdiff(colnames(meta_sub), c("specimenID", "tissue"))
  baseFormula <- "~1"

  mixed_vars <- intersect(variables,
                          c("sequencingBatch",  "rnaBatch", "libraryBatch", "cohort"))

  fixed_vars <- setdiff(variables, mixed_vars)

  if (length(mixed_vars) > 0) {
    mixed_vars <- paste0("(1 | ", mixed_vars, ")")
  }

  results <- mvIC::mvForwardStepwise(data_sub,
                                     baseFormula = baseFormula,
                                     data = meta_sub,
                                     variables = c(fixed_vars, mixed_vars))
  print(results$formula)

  saveRDS(results,
          file = file.path("data", "regression",
                           str_glue("{dataset}_{tissue}_mvIC_results.rds")))
}

# Mayo DLPFC: ~ rsem_uniquely_aligned_percent + RIN + percent_gc_content_fastq? + (1 | cohort)
# Mayo STG: ~ rsem_uniquely_aligned_percent + RIN + percent_gc_content_fastq? + (1 | cohort)
# Mayo CN: ~ rsem_alignable_percent + RIN + picard_PERCENT_DUPLICATION? + (1 | cohort)

# Could also remove multiqc / fastqc metrics, and add + (1 | rnaBatch)
