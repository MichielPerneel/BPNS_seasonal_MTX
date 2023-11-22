library(DESeq2)
# Specify data to be used
data_to_transform = 'data/kallisto/merged/tpm'


data <- read.csv(paste(data_to_transform, '.csv', sep = ''), # nolint
                    row.names = 1,
                    check.names = FALSE)

meta <- read.csv('samples.csv', sep = ';')

# Remove rows which are all zeros
data <- data[rowSums(data) > 0, ]

# Only perform below if counts and metadata have similar properties
# i.e. count columns == metadata rows
# If not, continue to next step
# dds <- DESeqDataSetFromMatrix(countData = round(data), colData = meta, design = ~ month)

# Add pseudocount if necessary
data_vst <- varianceStabilizingTransformation(as.matrix(round(data + 1)))
# data_vst <- varianceStabilizingTransformation(as.matrix(round(data)))

write.csv(data_vst, paste(data_to_transform, '_vst.csv', sep = ''), row.names = TRUE)
