# Load tidyverse package
library(tidyverse)

# Check if the correct number of command-line arguments is provided
if (length(commandArgs(trailingOnly = TRUE)) < 3) {
  cat("Usage: Rscript script_name.R pca_file eigenval_file location_file\n")
  quit(status = 1)
}

# Read command-line arguments
pca_file <- commandArgs(trailingOnly = TRUE)[1]
eigenval_file <- commandArgs(trailingOnly = TRUE)[2]
location_file <- commandArgs(trailingOnly = TRUE)[3]

# Read in PCA data
pca <- read_table2(pca_file, col_names = FALSE)
eigenval <- scan(eigenval_file)

# Print or plot the eigenvalues
print(eigenval)
plot(eigenval)

# Cleaning up the data
# Remove nuisance column
pca <- pca[,-1]
# Set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# Read location table
location_table <- read.table(location_file, header = TRUE, sep = "\t")

# Merge the PCA data with the location table based on individual IDs
pca <- merge(pca, location_table[, c("Sample", "Country")], by.x = "ind", by.y = "Sample", all.x = TRUE)

# Remake data.frame as a tibble
pca <- as_tibble(pca)
print(pca)

# Identify the numeric columns (excluding the first column 'ind' and the last column 'Country')
#numeric_columns <- setdiff(names(pca), c("ind", "Country"))

# Convert columns to numeric values
#pca[, numeric_columns] <- lapply(pca[, numeric_columns], function(x) as.numeric(as.character(x)))

#print(pca)


# Convert to percentage variance explained
pve <- data.frame(PC = 2:11, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = Country)) + geom_point(size = 3)
b <- b + scale_colour_discrete()
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

# Save the percentage variance explained plot
ggsave("percentage_variance_explained_plot.png", plot = a, width = 7, height = 5, units = "in")

# Save the PCA plot
ggsave("pca_plot.png", plot = b, width = 7, height = 5, units = "in")


