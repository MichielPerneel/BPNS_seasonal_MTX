library(ggplot2)
library(mgcv)
library(gridExtra)
library(dplyr)
library(grid)
library(tidygam)

# Load and preprocess data
data <- read.csv("data/analysis/WGCNA/pathway_tax_TPM_sum.csv")
total_pathway_data <- read.csv("data/analysis/WGCNA/pathway_TPM_sum.csv")

meta <- read.csv("samples.csv", sep = ";", header = TRUE, stringsAsFactors = FALSE)
meta$date <- as.Date(meta$date, format = "%d/%m/%Y")
meta$julian <- as.numeric(format(meta$date, "%j"))

# Merge core-Noctilucales and Dinophyceae into Dinophyceae
data$Taxogroup2_UniEuk[data$Taxogroup2_UniEuk == 'core-Noctilucales'] <- 'Dinophyceae'
# Sum the TPM values for core-Noctilucales and Dinophyceae
data <- data %>%
  group_by(Sample, Pathway_ID, Pathway_Name, Pathway_Class, Taxogroup2_UniEuk) %>%
  summarise(TPM = sum(TPM), .groups = 'drop')

# Add metadata
data <- data %>% left_join(meta, by = c("Sample" = "sample"))
total_pathway_data <- total_pathway_data %>% left_join(meta, by = c("Sample" = "sample"))

# Add year_month column
data$year_month <- format(data$date, "%Y-%m")
total_pathway_data$year_month <- format(total_pathway_data$date, "%Y-%m")

# Treat taxonomic groups as factors
data$Taxogroup2_UniEuk <- factor(data$Taxogroup2_UniEuk)

# Filter for selected taxonomic groups
data_sub <- data[data$Taxogroup2_UniEuk %in% c('Diatomeae', 'Dinophyceae', 'Prymnesiophyceae'),]

# Define color for each taxonomic group
tax_colors <- c("Diatomeae"="#E69F00", "Dinophyceae"="#56B4E9", "Prymnesiophyceae"="#009E73")

# Function to create plots for each pathway
create_plots <- function(pathway) {
  pathway_data <- subset(data_sub, Pathway_ID == pathway)
  pathway_total_data <- subset(total_pathway_data, Pathway_ID == pathway)
  pathway_annotated <- subset(data, Pathway_ID == pathway)

  if (length(unique(pathway_data$Taxogroup2_UniEuk)) > 1) {
    # GAM Plot
    gam_plot <- ggplot(pathway_data) +
      geom_point(aes(x=date, y=TPM, color=Taxogroup2_UniEuk)) +
      geom_smooth(aes(x=date, y=TPM, color=Taxogroup2_UniEuk), method="gam", formula=y ~ s(x)) +
      #geom_smooth(aes(x=date, y=TPM, color=Taxogroup2_UniEuk)) +
      coord_cartesian(xlim = as.Date(c("2020-07-01", "2021-08-01"))) +
      scale_color_manual(values = tax_colors) +
      scale_fill_manual(values = tax_colors) +
      labs(x="Date", y="TPM", color="Taxogroup2_UniEuk") +
      theme_minimal() +
      theme(legend.position="bottom") +
      xlim(as.Date("2020-07-01"), as.Date("2021-08-01"))

    # Marginal Plot
    first_sampling_dates <- total_pathway_data %>%
      group_by(year_month) %>%
      summarize(first_sampling_date = min(date))

    marginal_data <- merge(pathway_total_data, first_sampling_dates, by.x = "year_month", by.y = "year_month", all.x = TRUE)

    monthly_sum <- marginal_data %>%
      group_by(first_sampling_date) %>%
      summarise(monthly_TPM_sum = sum(TPM))

    # Same for all taxonomically annotated data
    marginal_data_annotated <- merge(pathway_annotated, first_sampling_dates, by.x = "year_month", by.y = "year_month", all.x = TRUE)

    monthly_sum_annotated <- marginal_data_annotated %>%
      group_by(first_sampling_date) %>%
      summarise(monthly_TPM_sum_annotated = sum(TPM))  # Check if you need TPM or another column here

    marginal_plot <- ggplot() +
      geom_bar(data=monthly_sum, aes(x=first_sampling_date, y=monthly_TPM_sum), stat="identity", fill="darkgrey") +
      geom_bar(data=monthly_sum_annotated, aes(x=first_sampling_date, y=monthly_TPM_sum_annotated), stat="identity", fill="#55af2b", alpha=0.5) +
      labs(x="", y="TPM") +
      theme_minimal() +
      ggtitle(unique(pathway_data$Pathway_Name))
    # Combine the two plots
    combined_plot <- grid.arrange(marginal_plot, gam_plot, ncol=1, heights=c(1, 1.5))
    # Save the combined plot
    ggsave(paste0("figures/WGCNA/pathways_GAM/gam/", pathway, "_plot.pdf"), combined_plot, width=6, height=4.5 )

  } else {
    warning(paste("Pathway", pathway, "has only one taxonomic group at the selected taxonomic level. No GAM plot will be created."))
  } 
}

# Apply the function to each unique pathway
lapply(unique(data$Pathway_ID), create_plots)
# Test
#create_plots("map00710")