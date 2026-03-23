# plot_tin_scores.R

args <- commandArgs(trailingOnly=TRUE)

# Read the input file path from command line arguments
input_file <- args[1]
output_file <- args[2]

# Load necessary packages
library(ggplot2)
library(reshape2)

# Load the data
tin_scores <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Reshape data from wide to long format
long_tin_scores <- melt(tin_scores, id.vars = "transcript")

# Rename columns
colnames(long_tin_scores) <- c("transcript", "sample", "TIN")

# Generate boxplots
p <- ggplot(long_tin_scores, aes(x = sample, y = TIN)) +
  geom_boxplot(fatten = 2, outlier.shape = NA, fill = "grey", color = "red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Sample", y = "TIN score", title = "TIN scores across samples")

ggsave(output_file, plot=p, width = 16, height = 12, units = "cm")