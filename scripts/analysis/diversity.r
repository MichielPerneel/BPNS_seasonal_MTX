library(dplyr)
library(vegan)

data_path <- "data/analysis/dinoflagellate_Name_to_Use_sample.csv"
output_path <- "data/analysis/dinoflagellate_biodiversity_estimates.csv"

data <- read.table(data_path,
                    header = TRUE, check.names = FALSE, sep = ",")
rownames(data) <- data[, 1]
data[, 1] <- NULL

# Remove all rows in which all values are 0
data <- data[rowSums(data) > 0, ]

# Transpose the data if necessary, should be genera x samples
# data <- t(data)

# normalize matrix: divide each row or column by the total of that row or column
data <- sweep(data, 1, rowSums(data), "/") # row
#data <- sweep(data, 2, colSums(data), "/") # column

shannon <- diversity(data, MARGIN = 2, index = "shannon")
simpson <- diversity(data, MARGIN = 2, index = "simpson")
beta_dist <- vegdist(data, method = "bray")
species <- specnumber(data, MARGIN = 2)
evenness <- shannon / log(species)

out <- data.frame(shannon, simpson, species, evenness)

write.csv(out, output_path, row.names = TRUE)
