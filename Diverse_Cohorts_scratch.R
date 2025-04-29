library(Matrix)
library(matrixStats)

source("helper_functions.R")

metadata <- download_metadata()

dataGen_fill <- viridis::turbo(3, begin = 0.25, end = 0.7)
# Same as above but minus #20202000 to blue and orange and minus #40404000 to
# green to darken
dataGen_color <- c("#089BCCFF", "#4FBF09FF", "#DE8612FF")

ggplot(metadata, aes(x = dataGenerationSite, y = RIN, fill = dataGenerationSite)) +
  geom_violin() + geom_boxplot(position = "dodge", width = 0.25) +
  scale_fill_manual(values = dataGen_fill) +
  theme_bw()

ggplot(metadata, aes(x = dataGenerationSite, y = DV200, fill = dataGenerationSite)) +
  geom_violin() + geom_boxplot(position = "dodge", width = 0.15) +
  scale_fill_manual(values = dataGen_fill) +
  theme_bw()

ggplot(metadata, aes(x = RIN, y = DV200, color = dataGenerationSite)) +
  geom_jitter(size = 0.3) +
  scale_color_manual(values = dataGen_color) +
  theme_bw()

columbia <- download_rsem("syn64289221")
mssm <- download_rsem("syn64176431")
mayo <- download_rsem("syn64176419")
rush <- download_rsem("syn64176441")

# Remove the 7 duplicate Rush samples
rush <- rush[, colnames(rush) %in% metadata$specimenID]

all(colnames(columbia) %in% metadata$specimenID)
all(colnames(mssm) %in% metadata$specimenID)
all(colnames(mayo) %in% metadata$specimenID)
all(colnames(rush) %in% metadata$specimenID)

columbia_log <- lognorm(columbia)
mssm_log <- lognorm(mssm)
mayo_log <- lognorm(mayo)
rush_log <- lognorm(rush)

variable_genes <- function(data, n_variable = 5000) {
  vv <- rowVars(data) |>
    sort(decreasing = TRUE)
  names(vv[1:n_variable])
}

# 2528 genes show up in all 4
vg_col <- variable_genes(columbia_log, 5000)
vg_mssm <- variable_genes(mssm_log, 5000)
vg_mayo <- variable_genes(mayo_log, 5000)
vg_rush <- variable_genes(rush_log, 5000)

pc_col <- prcomp(columbia_log[vg_col, ], scale. = TRUE)$rotation |>
  merge(metadata, by.x = "row.names", by.y = "specimenID")
pc_mssm <- prcomp(mssm_log[vg_mssm, ], scale. = TRUE)$rotation |>
  merge(metadata, by.x = "row.names", by.y = "specimenID")
pc_mayo <- prcomp(mayo_log[vg_mayo, ], scale. = TRUE)$rotation |>
  merge(metadata, by.x = "row.names", by.y = "specimenID")
pc_rush <- prcomp(rush_log[vg_rush, ], scale. = TRUE)$rotation |>
  merge(metadata, by.x = "row.names", by.y = "specimenID")


# Columbia ---------------------------------------------------------------------
# There are 2-3 distinct clusters, batch B08 is present in all of them so the
# separation isn't just by batch
ggplot(pc_col, aes(x = PC1, y = PC2, color = sequencingBatch)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "turbo")

ggplot(pc_col, aes(x = PC1, y = PC2, color = sequencingBatch)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  facet_wrap(~sequencingBatch)


# Clear separation by tissue: caudate nucleus is by itself, DLPFC and TP are
# mostly together
ggplot(pc_col, aes(x = PC1, y = PC2, color = tissue)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "viridis")


# All sample exchange samples are in the top cluster, but there are only 5 of them
ggplot(pc_col, aes(x = PC1, y = PC2, color = isSampleExchange)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  facet_wrap(~isSampleExchange)

# Not shown, no obvious separation by RIN, rnaBatch, libraryBatch, totalReads,
# cohort, sex, race, isHispanic, ageDeath, PMI, apoeGenotype,
# amyThal, Braak, amyCerad, amyA, amyAny, bScore, mayoDx, reag, ADoutcome


# MSSM -------------------------------------------------------------------------
# There's mostly one cluster with a long tail. Batches that aren't "B01" are
# mostly in the tail.

ggplot(pc_mssm, aes(x = PC1, y = PC2, color = sequencingBatch)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "turbo")

# Minor separation by tissue: caudate nucleus is mostly in the tail while
# DLPFC and STG are together
ggplot(pc_mssm, aes(x = PC1, y = PC2, color = sequencingBatch)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  facet_wrap(~sequencingBatch)

ggplot(pc_mssm, aes(x = PC1, y = PC2, color = tissue)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "viridis") +
  facet_wrap(~tissue)

# More "White" samples in the tail than in the main body of the cluster
ggplot(pc_mssm, aes(x = PC1, y = PC2, color = race)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  facet_wrap(~race)

# Samples with higher RIN closer to main body of cluster, really low RIN samples
# in the tail
ggplot(pc_mssm, aes(x = PC1, y = PC2, color = RINbin)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "viridis") +
  facet_wrap(~RINbin)

# Not shown: no obvious split by rnaBatch, libraryBatch, sex, isHispanic,
# ageDeath, PMI, amyCerad, Braak, amyAny, bScore, reag, ADoutcome, apoeGenotype


# Mayo -------------------------------------------------------------------------
# There are 2 clusters but no obvious separation by batch

ggplot(pc_mayo, aes(x = PC1, y = PC2, color = sequencingBatch)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "turbo")

ggplot(pc_mayo, aes(x = PC1, y = PC2, color = sequencingBatch)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  facet_wrap(~sequencingBatch)

# Clear separation by tissue: caudate nucleus is by itself while DLPFC and STG
# are together
ggplot(pc_mayo, aes(x = PC1, y = PC2, color = tissue)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "viridis")

# The white samples are mostly in the top cluster
ggplot(pc_mayo, aes(x = PC1, y = PC2, color = race)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "viridis") +
  facet_wrap(~race)

# Some bias for isHispanic with more FALSE values in the top cluster -- but this
# is correlated with race = White too
ggplot(pc_mayo, aes(x = PC1, y = PC2, color = isHispanic)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "viridis") +
  facet_wrap(~isHispanic)

# There might be a little bias toward contribution group (Emory vs Mayo), with
# Emory samples mostly in the top cluster but Mayo spread well between both
# clusters
ggplot(pc_mayo, aes(x = PC1, y = PC2, color = dataContributionGroup)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "viridis") +
  facet_wrap(~dataContributionGroup)

# *Possible* bias toward amyCerad ("None" are mostly in the top cluster) but the
# number of samples is probably too small to tell. Same with amyAny.
ggplot(pc_mayo, aes(x = PC1, y = PC2, color = amyCerad)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "viridis") +
  facet_wrap(~amyCerad)

# Some bias toward mayoDx and ADoutcome, with "Control" samples mostly in the
# top cluster. However this is a small number of samples. Same with reag ("No
# AD" and "Low Likelihood" in the top cluster).
ggplot(pc_mayo, aes(x = PC1, y = PC2, color = mayoDx)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "viridis") + facet_wrap(~mayoDx)

# Not shown: no obvious separation by RIN, rnaBatch, libraryBatch, sex, cohort
# (except for cohort == Emory), ageDeath, PMI, apoeGenotype, amyThal, Braak,
# amyA, bScore


# Rush -------------------------------------------------------------------------
# There are 2 distinct clusters. Batch B50 is only in the bottom cluster but
# the other two batches are pretty evenly in both.
ggplot(pc_rush, aes(x = PC1, y = PC2, color = sequencingBatch)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "turbo")

ggplot(pc_rush, aes(x = PC1, y = PC2, color = sequencingBatch)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  facet_wrap(~sequencingBatch)

# Clear separation by tissue: caudate nucleus is by itself while DLPFC and STG
# are together
ggplot(pc_rush, aes(x = PC1, y = PC2, color = tissue)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "viridis")

# All sample exchange samples are in the lower cluster, and most of them
# are NOT in batch B50
ggplot(pc_rush, aes(x = PC1, y = PC2, color = isSampleExchange)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  facet_wrap(~isSampleExchange)

# All "other" samples are in the bottom cluster and none of them are in B50.
# But this is a really small number of samples (16). No obvious split in other
# races.
ggplot(pc_rush, aes(x = PC1, y = PC2, color = race)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  facet_wrap(~race)

# Not shown: no obvious separation by RIN, rnaBatch (except for those in
# sequencingBatch == B50), libraryBatch (same), dataContributionGroup (except
# sample swaps), cohort, sex, isHispanic, ageDeath, PMI, apoeGenotype, amyThal,
# amyCerad, Braak, amyA, amyAny, bScore, reag, ADoutcome


# All data sets together -------------------------------------------------------

all_data <- cbind(columbia_log, mssm_log, mayo_log, rush_log)
metadata <- subset(metadata, specimenID %in% colnames(all_data))

table(metadata$tissue)

## Caudate nucleus

meta_caudate <- subset(metadata, tissue == "caudate nucleus")
data_caudate <- all_data[, meta_caudate$specimenID]

vg_caudate <- variable_genes(data_caudate, 5000)

pc_caudate <- prcomp(data_caudate[vg_caudate, ], scale. = TRUE)$rotation |>
  merge(meta_caudate, by.x = "row.names", by.y = "specimenID")

# Obvious batch effects by dataGenerationSite if you find variable genes on
# all data sets together
ggplot(pc_caudate, aes(x = PC1, y = PC2, color = dataGenerationSite)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "viridis")

ggplot(pc_caudate, aes(x = PC1, y = PC2, color = dataGenerationSite)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "viridis") +
  facet_wrap(~dataGenerationSite)

# Same story with dataContributionGroup
ggplot(pc_caudate, aes(x = PC1, y = PC2, color = dataContributionGroup)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "viridis")

ggplot(pc_caudate, aes(x = PC1, y = PC2, color = dataContributionGroup)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "viridis") +
  facet_wrap(~dataContributionGroup)

# And cohort
ggplot(pc_caudate, aes(x = PC1, y = PC2, color = cohort)) +
  geom_point() +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  facet_wrap(~cohort)


#rm(columbia_log, mayo_log, mssm_log, rush_log, all_data)
#gc()

# Use raw counts for GLM
all_data <- cbind(columbia, mssm, mayo, rush)
data_caudate <- all_data[, meta_caudate$specimenID]

# Some variables of interest
form <- paste(
  "~ amyThal + amyCerad + Braak + ageBin + pmiNumeric + race +",
  "isHispanic + RIN + (1|individualID) + (sequencingBatch | dataGenerationSite)"
) |>
  as.formula()

n_cores <- parallel::detectCores()-1

data_v <- variancePartition::voomWithDreamWeights(
  counts = as.matrix(data_caudate),
  formula = form,
  data = meta_caudate,
  plot = TRUE,
  save.plot = FALSE,
  BPPARAM = BiocParallel::SnowParam(n_cores)
)

fit_dream <- variancePartition::dream(
  exprObj = voom_counts,
  formula = form,
  data = meta_caudate,
  computeResiduals = TRUE,
  BPPARAM = BiocParallel::SnowParam(n_cores))

saveRDS(data_v, "caudate_voom.rds")
