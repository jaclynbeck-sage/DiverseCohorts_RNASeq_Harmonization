install.packages(c("BiocManager", "stringr", "dplyr", "ggplot2", "ggforce",
                   "patchwork", "viridis", "fastqcr"))

BiocManager::install(c("edgeR", "cqn"))

install.packages("synapser", repos = c("http://ran.synapse.org", "https://cloud.r-project.org"))
