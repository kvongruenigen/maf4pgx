################################################################################
#                         EDA of acquired data                                 #
################################################################################
#___________________________________________####################################
# Notes ------------------------------------------------------------------------

#'@Tasks:
## Describe the data and analyze what is possible to with the data
## Double Hits

#'@To-Dos:
# - Affected genes
# - Gene length normalization
# - Distinguish CNV DEL from SNV DEL
# - How to handle list information in R?
#   - snvs$clinical_interpretations: Extract list info
#   - snvs$consequence: Extract list info
# - Age per project
# - Samples per organ site
# - Mutations per chromosome
# - Extract clinical_information values
# - Find empty project in metadata (it is in database)
#   - biosample_id: pgxbs-kftvhnil
#   - Has only CNVs

#'@Done:
# - Data curation
# - Metadata summary -> Metadata overview



#'@3 data sets:
# - metadata (individuals + biosamples)
# - cnv variants
# - snv variants


#'@Tips

#'*Add mean and median to plot*'
# Add mean:
# geom_vline(aes(xintercept = mean(cnvfraction)),
#            color = "red",
#            linetype = "dashed",
#            size = 1) +
# annotate("text", x = mean(cnvs$cnvfraction),
#          y = 10,
#          label = paste("Mean =", round(mean(cnvs$cnvfraction), 2)),
#          color = "red") +
# 
# # Add median:
# geom_vline(aes(xintercept = median(cnvfraction)),
#            color = "blue",
#            linetype = "dashed",
#            size = 1) +
# annotate("text", x = median(cnvs$cnvfraction),
#          y = 15,
#          label = paste("Median =", round(median(cnvs$cnvfraction), 2)),
#          color = "blue") +
#

#'[Blue]

#'@Purple    

#'*Orange*

  

# Prepare ----------------------------------------------------------------------
rm(list = ls())
setwd("/Users/kayvongrunigen/Projects/maf4pgx/data")

library(tidyverse)
library(ggplot2)
library(ggmap)
library(ggpubr)
library(ggfortify)
library(GGally)
library(lubridate)
library(scales)
library(svglite)
library(tableone)
library(survival)
library(gt)
library(kableExtra)

# Themes -----------------------------------------------------------------------

# Custom color palette
my_color_palette <- c("#4C72B0", "#55A868", "#C44E52", "#8172B3", "#CCB974")

theme_basic <-   
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
        axis.text.x = element_text(angle = 0, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14))

theme_tight_bars <- 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14)
  )

theme_poster <-   theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 28),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 26),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 26)
  )

theme_alternative <- theme_classic +
  theme(
    text = element_text(family = "Times New Roman", size = 12),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 14),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",  # Adjust legend position as needed
    panel.grid.major = element_line(colour = "lightgray", linetype = "dashed", size = 0.2),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    plot.margin = margin(1, 1, 1, 1, "cm")  # Adjust plot margins as needed
  )


my_bar_chart_theme <- function() {
  theme_minimal() +
    theme(
      text = element_text(family = "Arial", size = 12, color = "black"),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 14),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "lightgray"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 12),
      plot.background = element_rect(fill = "white"),
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold", color = "black"),
      legend.background = element_rect(fill = "gray95"),
      legend.box.background = element_rect(color = "black"),
      legend.key = element_rect(color = "white"),
      legend.key.size = unit(1.2, "lines"),
      legend.margin = margin(5, 5, 5, 5),
      legend.spacing.x = unit(0.2, "lines"),
      legend.spacing.y = unit(0.2, "lines")
    )
}


# Functions --------------------------------------------------------------------
plot_top_levels <- function(data, variable, top_n = 10, include_others = TRUE, decreasing = TRUE) {
  # Convert the variable to a factor and sort levels by count in descending order
  data[[variable]] <- factor(
    data[[variable]],
    levels = names(sort(table(data[[variable]]), decreasing = decreasing))
  )
  
  # Create a table of counts and sort by count in descending order
  counts <- table(data[[variable]])
  sorted_levels <- names(sort(counts, decreasing = decreasing))
  
  # Identify the top levels
  top_levels <- sorted_levels[1:min(top_n, length(sorted_levels))]
  
  # Create a new factor variable, combining top levels and "Others"
  data$parsed_variable <- factor(
    data[[variable]],
    levels = c(top_levels, "Others")
  )
  
  # Assign the "Others" label to levels not in the top specified number
  data$parsed_variable[!(data$parsed_variable %in% top_levels)] <- "Others"
  
  # Filter the data to include only the top levels
  data_filtered <- data[data[[variable]] %in% top_levels, ]
  
  # Create a new factor variable for the original variable with levels as top levels
  data_filtered[[variable]] <- factor(
    data_filtered[[variable]],
    levels = top_levels
  )
  if (include_others){
    # Plot the ordered bar chart
    plot <- ggplot(data, aes(parsed_variable)) +
      geom_bar(fill = "lightblue", col = "black") +
      labs(title = paste("Top", top_n, "frequent", variable),
           x = "", y = "Count") +
      theme_basic
    
    return(plot)
  } else {
    # Plot the ordered bar chart
    plot <- ggplot(data_filtered, aes(parsed_variable)) +
      geom_bar(fill = "lightblue", col = "black") +
      labs(title = paste("Top", top_n, "frequent", variable),
           x = "", y = "Count") +
      theme_basic
    
    return(plot)
  }
  
}


# Data -------------------------------------------------------------------------
## Metadata --------------------------------------------------------------------

# Load data
metadata <- read_csv("progenetix/progenetix_tcga_individuals_biosamples_combined_data.csv")

# Order columns alphabetically for better overview
metadata <- metadata[, order(names(metadata))]

# Curate data
metadata <- metadata %>%
  # Convert multiple columns to factors
  mutate_at(vars(stage, substage, histological_diagnosis, icdo_morphology,
                 icdo_topography, sample_origin, tumor_type, project,
                 vital_status, sex, ethnicity),
            as.factor) %>%
  
  # Drop disease column
  # metadata$disease == metadata$histological_diagnosis
  # Expect Ductal carcinoma: histological_diagnosis is more precise
  select(-disease) %>%
  
  # Extract list: callset_ids
  unnest(callset_ids) %>%
  mutate(callset_ids = str_replace_all(callset_ids, "[\\[\\]']", ""))

## CNVs ------------------------------------------------------------------------

# Load data
cnvs <- read_csv("progenetix/progenetix_tcga_cnv_variants_data_with_genes.csv")

# Order columns alphabetically for better overview
cnvs <- cnvs[, order(names(cnvs))]

# Curate data
cnvs <- cnvs %>%
  # Convert multiple columns to factors
  mutate_at(vars(chromosome, cnv_state, relative_copy_class, variant_state),
            as.factor) %>%

  # Drop redundant column relative_copy_class & cnv_state
  # Relative_copy_class == cnv_state == variant_state
  # variant_state best phrasing ('copy number gain/loss')
  select(-c(relative_copy_class, cnv_state)) %>%

  # Extract list: affected_genes
  unnest(affected_genes) %>%
  mutate(affected_genes = str_replace_all(affected_genes, "[\\[\\]']", ""))

# Combine info from other data
## Add project from metadata to snvs over biosample_id
projects <- metadata %>%
  select(biosample_id, project)

cnvs <- cnvs %>%
  left_join(projects, by = "biosample_id")
rm(projects)


## SNVs ------------------------------------------------------------------------

# Load data
snvs <- read_csv("maf_analysis_data.csv")


# Set all columns to lower case
names(snvs) <- tolower(names(snvs))

# Separate allele frequency data frame
afs <- snvs %>%
  select(contains("af"))


# Curate data
snvs <- snvs %>%
  # Remove allele frequency columns
  select(-contains("af")) %>%
  
  # Rename column Hugo_Symbol to gene
  rename(gene = "hugo_symbol") %>%
  
  # Convert multiple columns to factors
  mutate_at(vars(biotype, canonical, center, chromosome, clin_sig, gene, hotspot,
                 impact, mirna, one_consequence, variant_classification,
                 variant_type),
            as.factor) %>%
  
  # Extract list info
  ## Aminoacid changes
  unnest(aminoacid_changes) %>%
  mutate(aminoacid_changes = str_replace_all(aminoacid_changes, "[\\[\\]']", "")) %>%
  ## Consequences
  unnest(consequence) %>%
  mutate(consequence = str_replace_all(consequence,"[\\[\\]']", ""))

# Order columns alphabetically for better overview
snvs <- snvs[, order(names(snvs))]

# ## Clinical interpretations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# library(jsonlite)
# 
# # Replace ' with "
# snvs$clinical_interpretations <- gsub("'", '"', snvs$clinical_interpretations)
# snvs$clinical_interpretations <- gsub('O"D', "O'D", snvs$clinical_interpretations)
# snvs$clinical_interpretations <- gsub('5"', "5'", snvs$clinical_interpretations)
# snvs$clinical_interpretations <- gsub(" None", ' "None"', snvs$clinical_interpretations)
# 
# snvs <- snvs %>%
#   mutate(clinical_interpretations = map(
#     clinical_interpretations,
#     ~if (!is.na(.)) jsonlite::fromJSON(.) else list()
#     ))
# 
# # Extract clinical_relevance values
# relevance_values <- snvs %>%
#   pull(clinical_interpretations) %>%
#   map(pluck, "clinical_relevance")


# Combine info from other data
## Add project from metadata to snvs over biosample_id
projects <- metadata %>%
  select(biosample_id, project)

snvs <- snvs %>%
  left_join(projects, by = "biosample_id")
rm(projects)

## Add cnv_fraction to snvs
cnv_fraction <- cnvs %>%
  select(biosample_id, cnvfraction) %>%
  distinct(biosample_id, .keep_all = TRUE)

snvs <- snvs %>%
  left_join(cnv_fraction, by = "biosample_id")
rm(cnv_fraction)

# Remove underscores from levels for all factor variables and first letter uppercase
for (i in 1:ncol(snvs)) {
  if (is.factor(snvs[[i]])) {
    levels(snvs[[i]]) <- gsub("_", " ", levels(snvs[[i]]))
  }
}

# Assign NO to the NA in canonical
if (!("NO" %in% levels(snvs$canonical))) {
  levels(snvs$canonical) <- c(levels(snvs$canonical), "NO")
}
snvs$canonical[is.na(snvs$canonical)] <- "NO"

snvs %>%
  filter(is.na(project))

#___________________________________________####################################
# Results ----------------------------------------------------------------------
#___________________________________________####################################
####################
plots_cnv_complete
plots_snvs_complete
####################
# Metadata ---------------------------------------------------------------------
# Sorted projects
plot_projects


## Population  -----------------------------------------------------------------
# Dashboard
individuals_summary

## Samples ---------------------------------------------------------------------
sample_summary


## Metadata overview -----------------------------------------------------------

# ggarrange(individuals_summary, sample_summary, ncol = 1)

# ggarrange(
#   plot_sex,
#   plot_stage,
#   plot_ethnicities,
#   plot_tumor_type,
#   plot_age_by_sex,
#   plot_vital_status,
#   plot_followup_time,
#   plot_collection_moment,
#   plot_frequent_origins,
#   plot_frequent_histologicals,
#   plot_frequent_topographies,
#   plot_frequent_morphologies,
#   ncol = 4, nrow = 3
#   )

metadata_summary <-
ggarrange(
  ggarrange(
    plot_sex,
    plot_stage,
    plot_ethnicities,
    plot_tumor_type,
    plot_followup_time,
    ncol = 5, nrow = 1),
  
  ggarrange(
    plot_age_by_sex,
    plot_vital_status,
    ncol = 2, nrow = 1),
  
  ggarrange(
    plot_collection_moment,
    plot_frequent_origins,
    plot_frequent_histologicals,
    plot_frequent_topographies,
    plot_frequent_morphologies,
    ncol = 5, nrow = 1),
  
  ncol = 1, nrow = 3,
  heights = c(1, 2, 1)
)

ggsave("individuals_sample_summary.png",
       path = "../analysis/plots",
       plot = metadata_summary,
       width = 20, height = 25)

plot_age_per_project


## CNVS -------------------------------------------------------------------------

plots_cnv_complete

ggsave("cnv_summary.png",
       path = "../analysis/plots",
       plot = cnv_summary,
       width = 20, height = 25)



## SNVS ------------------------------------------------------------------------
# Sumamry snvs
plots_snvs_complete

# Variant types
plot_variant_types

# Variant classification
plot_variant_classification

# Impact of Variants
plot_impact

# Variants per sample and project
plot_variations_per_sample
plot_variations_per_sample_different_distr

# Top Consequences
plot_top3_consequenes

# Top 10 mutated genes
top10_mut_genes


ggarrange(
  ggarrange(
    plot_variant_types,
    plot_variant_classification,
    plot_impact,
    top10_mut_genes,
    ncol = 4, nrow = 1),
  plot_variants_per_sample,
  ncol = 1, nrow = 2)


## Together --------------------------------------------------------------------

# Variants vs CNV fraction
plot_snv_cnv_fraction

# Number of genes affected by CNVs and SNVs
plot_density_of_overlaps

#___________________________________________####################################
# Experiments ------------------------------------------------------------------
#___________________________________________####################################

# Oncoplot and Manhattan-style plot
library(ComplexHeatmap)
library(circlize)
library(maftools)

# Assuming cnvs and snvs have columns 'start', 'end', and 'type' (CNV or SNV)
cnv_positions <- cnvs %>% select(start, end, chromosome)
snv_positions <- snvs %>% select(start, end, chromosome)

# Create a data frame for oncoplot
oncoplot_data <- rbind(cnv_positions, snv_positions)
oncoplot_data$type <- rep(c("CNV", "SNV"), each = nrow(oncoplot_data) / 2)

# Convert the data to a format suitable for ComplexHeatmap
oncoplot_matrix <- Heatmap(as.matrix(oncoplot_data[, c("start", "end")]),
                           col = c("CNV" = "blue", "SNV" = "red"),
                           cluster_columns = FALSE,
                           show_column_names = FALSE,
                           name = "Type",
                           row_title = "Position")

# Draw the oncoplot using both ComplexHeatmap and circlize
draw(oncoplot_matrix)

# Create a Manhattan-style plot with chromosomal grouping
ggplot(oncoplot_data, aes(x = (start + end) / 2, y = significance, color = type)) +
  geom_point(size = 2) +
  facet_wrap(~chromosome, scales = "free_x", ncol = 1) +
  scale_color_manual(values = c("CNV" = "blue", "SNV" = "red")) +
  labs(title = "Manhattan-style Plot for CNVs and SNVs",
       x = "Genomic Position",
       y = "Significance") +
  theme_minimal()


ggplot(cnvs, aes(x = cnvfraction, y = variant_length)) +
  geom_point(aes(col = variant_state)) +
  labs(title = "CNV Fraction vs. Variation Length",
       x = "CNV Fraction",
       y = "Variation Length") +
  theme_basic

