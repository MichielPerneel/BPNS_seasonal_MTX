library(vegan)
library(mixOmics)
library(svglite)
library(pheatmap)
library(Hmisc)
library(RColorBrewer)

# Here we try to relate taxonomic information to functional diversity.

# ------- Variance partitioning -------
# First we tried variance partitioning, but results were off.
# The idea was to assess how much of the variation in the functional data,
# either in the form of eigengenes or in the form of functional identifier expression,
# could be attributed to the variation in the community composition.
# The varpart function from the mixOmics package was used for this, but the results were not very convincing.
# Best explanation we could think off why it didn"t work was that the functional data
# and the community data were not independent, but drawn from the same metatranscriptomic distribution (count data),
# which made comparing the two very difficult. 
# Also, the eigengene set was way smaller and so very likely to be fully explained by the larger community data.

# Then we abandond variance partitioning, and use PLS to look at the correlation between the community data and the eigengene expression.
# The goal was to try to infer which groups of taxa were correlated with the variation in the eigengene expression.

# ------- PLS -------
# Instead of variance partitioning, we used sPLS to look at the correlation between the community data and the eigengene expression.
# The goal was to try to infer which groups of taxa were correlated with the variation in the eigengene expression.
# See: http://mixomics.org/graphics/cim/ for a nice explanation of the PLS algorithm.

# Load data
eigengene <- read.csv("data/analysis/WGCNA/module_eigengenes.csv",
                        header = T, row.names = 1)
dim(eigengene)

# PhyloDB annotations
#taxonomy <- read.csv("data/analysis/sample_class_TPM_sum.csv",
#                        header = T, row.names = 1)
#taxonomy <- t(read.csv("data/analysis/species_sample_TPM_sum.csv",
#                        header = T, row.names = 1))

# Eukprot annotations
#taxonomy <- t(read.csv("data/analysis/eukprot_species_sample_TPM_sum.csv",
#                        header = T, row.names = 1, check.names = F))

taxonomy <- t(read.csv("data/analysis/eukprot_genus_sample_TPM_sum.csv",
                        header = T, row.names = 1, check.names = F))

#taxonomy <- t(read.csv("data/analysis/eukprot_order_sample_TPM_sum.csv",
#                        header = T, row.names = 1, check.names = F))
dim(taxonomy)

# Preprocess the data
## Check for missing values
any(is.na(eigengene))
any(is.na(taxonomy))

all(rownames(eigengene) == rownames(taxonomy))

## Delete low-abundant components
taxonomy <- taxonomy[,colSums(taxonomy) > 1]

### Plot distribution of column sums
plot(log(sort(colSums(taxonomy))), type = "h",
                            xlab = "Class",
                            ylab = "log(sum)")
### Add horizontal lines at points of interest
abline(h = log(50), col = "red")
abline(h = log(500), col = "red")
abline(h = log(1000), col = "red")
abline(h = log(5000), col = "red")
abline(h = log(50000), col = "red")

### Plot histogram of column sums
hist(log(colSums(taxonomy)),
        xlab = "log(sum)",
        ylab = "Frequency")
abline(v = log(50), col = "red")
abline(v = log(500), col = "red")
abline(v = log(1000), col = "red")
abline(v = log(5000), col = "red")
abline(v = log(50000), col = "red")

## Log transform taxonomy
taxonomy_log <- log(taxonomy + 1)
taxonomy_scaled <- scale(taxonomy)

# Calculate correlations and p-values
cor <- cor(taxonomy_log, eigengene, method = "pearson")

## Calculate p-values
pval <- matrix(NA, nrow = nrow(cor), ncol = ncol(cor))
for (i in 1:ncol(cor)) {
    for (j in 1:nrow(cor)) {
        pval[j, i] <- cor.test(
            as.numeric(taxonomy_log[, j]), as.numeric(eigengene[, i]))$p.value
        }
    }

# Adjust p-values for multiple testing
qval <- matrix(p.adjust(
    as.vector(pval), method = "BH"),
    nrow = nrow(pval),
    ncol = ncol(pval))

## Name rows and columns
rownames(qval) <- rownames(cor)
colnames(qval) <- colnames(cor)

# Order by abundance
taxonomy_abundance <- colSums(taxonomy)
taxonomy_ordered <- taxonomy_abundance[
    order(taxonomy_abundance, decreasing = TRUE)]

# Select the 30 most abundant taxonomic classes
top_30 <- head(taxonomy_ordered, 30)
cor <- cor[names(top_30), ]
qval <- qval[names(top_30), ]

# Check if column and row order is equal between cor and qval
all(rownames(cor) == rownames(qval))
all(colnames(cor) == colnames(qval))

# Save figure as svg (7 by 8 cm)
breaks_list <- seq(-0.75, 1, 0.01)

svglite("figures/biodiversity_functional_diversity/eukprot_genus_eigengene_correlation.svg", width = 7, height = 8)

pheatmap(cor, fontsize = 12,
        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaks_list)),
        breaks = breaks_list)
dev.off()

# Print which classes are significantly correlated with which eigengenes
for (i in 1:ncol(cor)) {
    print(colnames(cor)[i])
    print(rownames(cor)[which(qval[,i] < 0.05)])
}

## Heatmap should have asterisks for p-values < 0.05!
## Add asterisks to heatmap manually for the printed classes and eigengenes above

# PLS
## Preliminary analysis with PCA
pca_taxonomy <- pca(taxonomy_log,
                    ncomp = 10,
                    center = TRUE,
                    scale = TRUE)

plot(pca_taxonomy)

# PLS and tuning
pls_result <- pls(taxonomy_log,
                    eigengene,
                    mode = "regression")

# Save figure as svg
svglite("figures/biodiversity_functional_diversity/class_eigengene_PLS.svg")
par(mar = c(6, 8.5, 3, 3));
cim(pls.result)
dev.off()

perf_pls_result <- perf(pls_result,
                        validation = "Mfold",
                        folds = 10,
                        nrepeat = 2)

plot(perf_pls_result, criterion = "Q2.total")

# set range of test values for number of taxonomic categories to use from taxonomic data
list.keepX <- c(seq(1, 50, 5))

tune_spls_result <- tune.spls(taxonomy_log,
                            eigengene,
                            ncomp = 3,
                            test.keepX = list.keepX,
                            nrepeat = 1, folds = 10,
                            mode = "regression",
                            measure = "cor")
plot(tune_spls_result)

# extract optimal number of variables for X dataframe
optimal_keepX <- tune_spls_result$choice.keepX

optimal_ncomp <-  length(optimal_keepX) # extract optimal number of components

# use all tuned values from above
final_spls_result <- spls(taxonomy_log, eigengene,
                    ncomp = optimal_ncomp,
                    keepX = optimal_keepX,
                    mode = "regression")

perf_spls_result <- perf(final.spls.result,
                        validation = "Mfold",
                        folds = 10,
                        nrepeat = 2)
plot(perf_spls_result, criterion = "Q2.total")


plotVar(final_spls_result, cex = c(3,4), var.names = c(FALSE, TRUE))
# Save figure as svg
svglite("figures/biodiversity_functional_diversity/class_eigengene_PLS_tuned.svg", width = 10, height = 10)
par(mar = c(6, 8.5, 3, 3));
cim(final_spls_result)
dev.off()