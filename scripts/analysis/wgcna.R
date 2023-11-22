# WGCNA pipeline adapted from Natalie Cohen
# (https://github.com/cnatalie/METZYME/blob/master/WGCNA.R)

# Load libraries
library(WGCNA)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(flashClust)
library(DESeq2)
library(clusterProfiler)
library(svglite)
library(reshape2)
library(tidyr)
library(forcats)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

#Read in normalized, gene-annotated (modules or KOs) transcript counts
samples_genes_matrix <- read.csv("data/analysis/KEGG_ko_TPM_sum_transposed.csv",
                                    check.names = F)
                                    
# Set row names and column names
rownames(samples_genes_matrix) <- samples_genes_matrix[, 1]
samples_genes_matrix <- samples_genes_matrix[, -1]

# Set values below 1 to 0 (lower detection limit)
samples_genes_matrix[samples_genes_matrix < 1] <- 0

# Check the distribution of KO colsums
hist(log(colSums(samples_genes_matrix)), breaks = 100)
abline(v = 2, col = "red")

# Remove columns with colsums below 10
samples_genes_matrix <- samples_genes_matrix[, colSums(samples_genes_matrix) > 10]

# Transform to datExpr matrix
datExpr <- as.matrix(samples_genes_matrix)
#datExpr <- as.matrix(t(samples_genes_matrix))

# Check for outliers in the data
gsg <- goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK #if TRUE, no outlier genes

# log2 transform, add pseudocount
datExpr <- log2(datExpr + 1)

#Read in environmental data
data_env <- read.csv("data/environmental/samples_env.csv",
                    stringsAsFactors = F, sep = ";",
                    header = T, row.names = 1)

# The below check should return TRUE if datasets align correctly,
# if not the names are out of order
table(rownames(data_env) == rownames(datExpr))

# Get the metadata from env_data
metadata <- data_env[, c("month", "station")]

# Get the environmental parameters from env_data
env_params <- data_env[, c("NO3", "PO4", "Si", "SPM", "Temperature", "salinity")]

#-------------  1. Dendrogram and trait heatmap showing outliers ------------#
## Calculate sample network adjacency matrix
A <- adjacency(t(datExpr), type = "signed")
# type = signed => both positive and negative correlations are considered

k <- as.numeric(apply(A, 2, sum)) - 1
# k is a numeric vector that contains the sum of each column in A, subtracted by 1

Z.k <- scale(k)
# Z.k is a numeric vector that contains the scaled k values.
# Scale is a function that centers and scales the data.

thresholdZ.k <- -3
# Choose a threshold value (often -2.5)
# We took -3 since outliers samples had biological meaning (Phaeocystis bloom),
# and we wanted to keep them in the analysis.
outlierColor <- ifelse(Z.k < thresholdZ.k, "red", "black")
# Mark outliers in red

# Calculate the sample dendrogram
sampleTree <- flashClust(as.dist(1 - A), method = "average")

# Create a data frame with environmental parameters as colors for the heatmap
traitColors <- data.frame(numbers2colors(env_params, signed = FALSE))
dimnames(traitColors)[[2]] <- paste(names(env_params))

# Combine the outlier colors with the trait colors
datColors <- data.frame(outlierC = outlierColor, traitColors)

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,
                    groupLabels = names(datColors),
                    colors = datColors,
                    main = "Sample dendrogram and trait heatmap")

# Remove outlying samples from the expression matrix and trait matrix
samples_to_remove <- Z.k < thresholdZ.k | is.na(Z.k)
datExpr <- datExpr[!samples_to_remove, ]
env_params <- env_params[!samples_to_remove, ]
metadata <- metadata[!samples_to_remove, ]

#-------------  2. Network construction and module detection  ------------#
# Choose a set of soft-thresholding powers
powers <- c(seq(1, 20, by = 1), seq(20, 60, by = 5));

# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr,
                        powerVector = powers,
                        networkType = "signed",
                        verbose = 2) # we want smallest value, closest to 0.9
sizeGrWindow(3, 5)
par(mfrow = c(2, 1));
cex1 <- 0.9;
#Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Soft Threshold (power)",
    ylab = "Scale Free Topology Model Fit,signed R^2",
    type = "n",
    main = paste("Scale independence"));
text(sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    labels = powers, cex = cex1, col = "red");
abline(h = 0.80, col = "red")
abline(h = 0.90, col = "red")
#This line corresponds to using an R^2 cut-off of h

#Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)",
    ylab = "Mean Connectivity",
    type = "n",
    main = paste("Mean connectivity"))
text(sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    labels = powers,
    cex = cex1, col = "red")

# Choose soft thresholding power based on the plot above
softPower <- 14 # smallest value to plateau above 0.8
# Calculate the adjacency matrix, or re-use a previously calculated adjacency matrix
# adjacency <- attach("data/analysis/WGCNA/adjacency.RData")$adjacency
adjacency <- adjacency(datExpr,
                        power = softPower,
                        type = "signed")

# Store the adjacency matrix locally for future use, calculation takes a while
# save(adjacency, file = "data/analysis/WGCNA/adjacency.RData")
                        
# Translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

# Each branch corresponds to a KO id
# branches grouping together densely are co-expressed KOs
geneTree <- flashClust(as.dist(dissTOM), method = "average")
minModuleSize <- 70
dynamicMods <- cutreeDynamic(dendro = geneTree,
                            distM = dissTOM,
                            deepSplit = 4,
                            pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)
dynamicColors <- labels2colors(dynamicMods)

sizeGrWindow(3, 3)
plotDendroAndColors(geneTree,
                    dynamicColors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Merge modules whose expression profiles are very similar
MEList <- moduleEigengenes(datExpr,
                        colors = dynamicColors,
                        softPower = softPower)
MEs <- MEList$eigengenes

# Calculate dissimilarity of module eigenegenes
MEDiss <- 1 - cor(MEs)

# Cluster module eigengenes
METree <- flashClust(as.dist(MEDiss), method = "average")

sizeGrWindow(7, 6)
plot(METree,
    main = "Clustering of module eigengenes",
    xlab = "", sub = "")

MEDissThres <- 0.4
abline(h = MEDissThres, col = "red")

merge <- mergeCloseModules(datExpr,
                        dynamicColors,
                        cutHeight = MEDissThres,
                        verbose = 3)

mergedColors <- merge$colors
mergedMEs <- merge$newMEs

plotDendroAndColors(geneTree,
                    cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                        guideHang = 0.05)

moduleColors <- mergedColors
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder) - 1
MEs <- mergedMEs

#-------- 3. Relating modules to traits and finding important genes -------#
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
table(moduleColors)

# Recalculate module eigengenes with color labels
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
# Remove the grey module,
# The grey module contains all genes not assigned to any module
MEs0 <- subset(MEs0, select = -MEgrey)
MEs <- orderMEs(MEs0)
# Check correlation between module eigengenes and environmental parameters
moduleTraitCor <- cor(MEs, env_params, use = "p");
# Calculate p-values for the correlations
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples);

# Represent module trait correlations as a heatmap
sizeGrWindow(10, 6)
# Display correlations and their p-values
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) <- dim(moduleTraitCor)


# Create a heatmap plot
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
             xLabels = names(env_params),
             yLabels = names(MEs),
             ySymbols = names(MEs),
             colorLabels = FALSE,
             colors = blueWhiteRed(50),
             textMatrix = textMatrix,
             setStdMargins = FALSE,
             cex.text = 0.8,
             zlim = c(-1, 1),
             main = "Module-trait relationships")
# Save the plot
svglite("figures/WGCNA/module_trait_relationships.svg")
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
             xLabels = names(env_params),
             yLabels = names(MEs),
             ySymbols = names(MEs),
             colorLabels = FALSE,
             colors = blueWhiteRed(50),
             textMatrix = textMatrix,
             setStdMargins = FALSE,
             cex.text = 0.8,
             zlim = c(-1, 1),
             main = "Module-trait relationships")
dev.off()

# Save the module eigengene expression data as csv
write.csv(MEs, file = "data/analysis/WGCNA/module_eigengenes.csv", row.names = TRUE)

# Plot the expression of the eigengene of each module at a given month and station using boxplots
## First add the index as a column
MElong <- MEs %>% as.data.frame() %>% rownames_to_column("sample")
## Reshape the data to the long format, add station and month as factors
MElong <- melt(MElong,
            id_vars = "sample",
            variable.name = "module",
            value.name = "ME_expression")

## metadata index should be column named sample
metadata2 <- metadata %>% rownames_to_column("sample")
## Add station and month metadata
MElong <- MElong %>% 
    left_join(metadata2, by = "sample")

# Plot module eigengene expression at each station and month as a heatmap
## Order the months
MElong$month <- factor(MElong$month, 
                        levels = c("July_2020", "August_2020", "September_2020",
                            "November_2020", "December_2020", "January_2021",
                            "February_2021", "April_2021", "May_2021",
                            "June_2021", "July_2021"))
## Order the stations
MElong$station <- factor(MElong$station,
                        levels = c("700", "780", "130", "330", "120", "ZG02"))

# Plot
ggplot(MElong, aes(x = month, y = station, fill = ME_expression)) +
    geom_tile() +
    facet_wrap(. ~ module) +
    scale_fill_distiller(palette = "Spectral") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Month", y = "Station", fill = "Module eigengene expression") +
    theme(legend.position = "bottom", ) +
    # Add vertical legend in the bottom right corner
    guides(fill = guide_colorbar(barwidth = 10, barheight = 1, title.position = "top", title.hjust = 0.5),
            colour = guide_legend(direction = "vertical", nrow = 1, byrow = TRUE),
            position = "bottom")
# Save as svg
ggsave("figures/WGCNA/module_eigengene_expression_per_month_station.svg", width = 8, height = 6)

#----------------------- 4. Module content ----------------------#
# For a given module, we define a quantitative measure of the importance of a KEGG ID in the module,
# by correlating the module eigengene with the KEGG ID expression (module membership)

KEGG_ID_ModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"));

# per module, rank the most correlated KEGG IDs in decreasing order of correlation
# For these genes, construct a dataframe that contains the KEGG ID name, the module membership, and the description
# The resulting dataframe is then used to find the most important KEGG IDs in each module through Mann-Whitney U tests

for (module in names(MEs)) {
    print(module)
    # Select the column of KEGG IDs correlating with the module from the dataframe
    # And order them by decreasing correlation
    module_cor <- select(KEGG_ID_ModuleMembership, module)
    module_cor <- module_cor[order(-module_cor[, 1]), , drop = FALSE]
    # Add the index as a column
    module_cor <- module_cor %>% as.data.frame() %>% rownames_to_column("KEGG_ID")
    print(nrow(module_cor))
    # Save the dataframe as a text file
    write.table(module_cor, 
                file = paste0("data/analysis/WGCNA/information/", module, "_content.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    }

# Print all KEGG IDs to a text file
write.table(rownames(KEGG_ID_ModuleMembership),
            file = "data/analysis/WGCNA/information/KEGG_ID_list.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)