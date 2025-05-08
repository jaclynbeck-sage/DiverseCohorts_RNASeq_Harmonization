install.packages(c("BiocManager", "stringr", "dplyr", "ggplot2", "ggforce",
                   "patchwork", "viridis", "fastqcr"))

BiocManager::install(c("biomaRt", "edgeR", "GenomicRanges", "rtracklayer",
                       "BSgenome", "cqn"))
# rtracklayer requires install of libbz2-dev

install.packages("synapser", repos = c("http://ran.synapse.org", "https://cloud.r-project.org"))
