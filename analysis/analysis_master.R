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

# Colors:
# https://stackoverflow.com/questions/50321000/r-colors-many-distinctive-colors-that-are-still-pretty


## Plot legend seperately ------------------------------------------------------
# legend <- get_legend(
#   ggplot(plot.long.data,aes(y=value))+
#     geom_boxplot(outlier.size=0.15,aes(fill=method))+
#     scale_fill_manual(values=cols)+ylim(c(0,1))+ylab('F1')
# )


## EDA Recipe: -----------------------------------------------------------------
# # 1. Independence of Y: Plot Y against time/space variables. (Not really in this data)
# # 2. Distribution of Y: 
# #   - Plot Y histogram as density with overlaid distribution
# #   - Make QQ plots
# 
# # Histograms
# ggarrange(
#   ggplot(mtcars, aes(disp)) +
#     geom_histogram(aes(y= ..density..), fill= "white", col= "black") + #density instead of counts
#     stat_function(fun= dnorm, #function for normal distribution
#                   args= list(mean(mtcars$disp),# R needs mean
#                              sd(mtcars$disp)), # and SD
#                   n= 1e2,
#                   col= "blue") +
#     labs(title= "Density histogram",
#          subtitle= "Gaussian distribution overlaid") +
#     theme_bw() +
#     theme(plot.title= element_text(hjust= 0.5),
#           plot.subtitle= element_text(hjust= 0.5)),
#   ggplot(mtcars, aes(disp)) +
#     geom_histogram(aes(y= ..density..), fill= "white", col= "black",
#                    binwidth= 1) +
#     stat_function(fun= dpois, # function for Poisson distribution
#                   args= list(mean(mtcars$disp)),
#                   n= max(mtcars$disp) - min(mtcars$disp) + 1,
#                   fill= "blue", col= NA, alpha= .2, geom= "bar") +
#     labs(title= "Density histogram",
#          subtitle= "Poisson distribution overlaid") +
#     theme_bw() +
#     theme(plot.title= element_text(hjust= 0.5),
#           plot.subtitle= element_text(hjust= 0.5)),
#   ncol= 2)
# 
# 
# # QQplots (possible to draw for different distributions)
# ggarrange(
#   ggplot(metadata, aes(sample= age_at_diagnosis)) +
#     stat_qq(distribution= stats::qpois,
#             dparams= list(lambda= mean(metadata$age_at_diagnosis, na.rm = T))) +
#     stat_qq_line(distribution= stats::qpois,
#                  dparams= list(lambda= mean(metadata$age_at_diagnosis, na.rm = T)), col= "red") +
#     labs(x= "Theoretical quantiles (Poisson)", y= "Sample Quantiles") +
#     theme_bw(),
#   ggplot(metadata, aes(sample= age_at_diagnosis)) +
#     stat_qq() +
#     stat_qq_line(col= "red") +
#     labs(x= "Theoretical quantiles (Gaussian)", y= "Sample Quantiles") +
#     theme_bw(),
#   ncol= 2)
# 
# # 3. Homogeneity of Y: scatterplots and conditional boxplots (use facet_wrap)
# # Across the covariate
# ggplot(owls, aes(ArrivalTime, NCalls)) +
#   geom_point() +
#   facet_grid(FoodTreatment~ SexParent) +
#   theme_bw()
# 
# # For each of the categorical predictor combinations, e.g.
# ggplot(owls, aes(FoodTreatment, NCalls)) +
#   geom_boxplot(fill= "lightgrey") +
#   facet_wrap(~ SexParent) +
#   theme_bw()
# 
# # 4. Outliers in Y and X: Boxplots and Cleveland dotplots
# ggarrange(
#   ggplot(owls, aes(y= NCalls)) +
#     geom_boxplot(fill= "lightgrey") +
#     labs(x= "", y= "Number of Calls") +
#     theme_bw() +
#     theme(axis.text.x= element_text(colour= NA),
#           axis.ticks.x= element_blank()),
#   ggplot(owls, aes(x= 1:nrow(owls), y= NCalls)) +
#     geom_point(alpha= .2) +
#     labs(x= "Order of the data", y= "Number of Calls") +
#     theme_bw(),
#   ggplot(owls, aes(y= ArrivalTime)) +
#     geom_boxplot(fill= "lightgrey") +
#     labs(x= "", y= "Arrival time") +
#     theme_bw() +
#     theme(axis.text.x= element_text(colour= NA),
#           axis.ticks.x= element_blank()),
#   ggplot(owls, aes(x= 1:nrow(owls), y= ArrivalTime)) +
#     geom_point(alpha= .2) +
#     labs(x= "Order of the data", y= "Arrival time") +
#     theme_bw(),
#   ncol= 2, nrow= 2)
# 
# ## v) Relationships among Y and X 
# ##    (entails the "Collinearity X", "Relationships Y & X" and "Interactions"
# ##    steps in the Zuur et al 2010 protocol; scatterplots, coplots and boxplots)
# 
# # Scatterplot-matrix (both functions here are highly customizable)
# # Default
# ggpairs(metadata, columns = c("age_at_diagnosis", "followup_time_months", "days_to_death"))
# 
# 
# # Alternatively
# scatterplotMatrix(owls[c(5:8)], smooth= F)
# 
# # Coplot, non-parallel lines might suggest significant interactions
# # (note: same plot used to assess homogeneity of Y across the covariate!)
# ggplot(owls, aes(ArrivalTime, NCalls)) +
#   geom_point() +
#   geom_smooth(method= "lm") +
#   facet_grid(FoodTreatment~ SexParent) +
#   theme_bw()


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
library(pastecs)
library(car)
library(dplyr)
library(survival)
library(survminer)
library(mongolite)
library(e1071)
library(pals)

# Themes -----------------------------------------------------------------------

# Custom color palette
my_color_palette <- c("#4C72B0", "#55A868", "#C44E52", "#8172B3", "#CCB974",
                      "#64B5CD", "#8172B3", "#937860", "#DA8BC3", "#8C8C8C")

max_colors <- polychrome()
names(max_colors) <- NULL

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

theme_poster <-
  theme_classic() +
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

# Make a function to calculate the standard error SE
se <- function(x) {
  sd(x, na.rm = TRUE) / sqrt(length(x))
}

# Data -------------------------------------------------------------------------
## Genes -----------------------------------------------------------------------
# Connect to Progenetix database to get gene length
gene_collection <- mongo(
  collection = "genes",
  db = "progenetix",
  url = "mongodb://localhost",
  verbose = FALSE,
  options = ssl_options()
)

# Add gene length to gene_counts data frame
# gene_locus_length = gene length
# symbol = gene name

# Get data table with gene length and gene name
genes <- gene_collection$find(
  fields = '{"symbol": 1, "gene_locus_length": 1, "end": 1, "start": 1, "orientation": 1, "reference_name": 1, "_id": 0}'
)

# Rename columns
genes <- genes %>%
  rename(
    gene = symbol,
    gene_start = start,
    gene_end = end,
    gene_orientation = orientation,
    gene_length = gene_locus_length,
    chromosome = reference_name)

genes

cgc <- read_tsv("Cosmic_CancerGeneCensus_v99_GRCh38.tsv")

# Rename columns
cgc <- cgc %>%
  rename(gene = GENE_SYMBOL,
         full_gene_name = NAME,
         cgc_chromosome = CHROMOSOME,
         cgc_start = GENOME_START,
         cgc_end = GENOME_STOP,
         cytoband = CHR_BAND,
         somatic = SOMATIC,
         germline = GERMLINE,
         tumour_types = TUMOUR_TYPES_SOMATIC,
         germline_tumour_types = TUMOUR_TYPES_GERMLINE,
         transloction_partner = TRANSLOCATION_PARTNER,
         gene_type = MOLECULAR_GENETICS,
         cancer_role = ROLE_IN_CANCER,
         cancer_syndrome = CANCER_SYNDROME,
         mutation_types = MUTATION_TYPES,
         tier = TIER
  ) %>%
  select(-c(COSMIC_GENE_ID, TISSUE_TYPE, OTHER_GERMLINE_MUT, OTHER_SYNDROME,
            SYNONYMS))


# Merge
genes <- genes %>%
  left_join(cgc, by = "gene")

# Show differences
genes %>%
  filter(!is.na(cgc_start)) %>%
  select(gene_end, cgc_end, gene_start, cgc_start,
         chromosome, cgc_chromosome) %>%
  # Calculate difference between start and end
  mutate(end_diff = gene_end - cgc_end,
         start_diff = gene_start - cgc_start,
         chrom = chromosome == cgc_chromosome)

genes <- genes %>%
  select(-c(cgc_start, cgc_end, cgc_chromosome))

head(genes)

cgc_genes <- genes %>%
  filter(!is.na(tier))

head(cgc_genes)



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

unique(metadata$icdo_morphology)
unique(metadata$icdo_topography)
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


# Merge with cgc cancer role
cgc$gene <- as.factor(cgc$gene)
snvs <- snvs %>%
  left_join(select(cgc, gene, cancer_role), by = "gene")



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
    #plot_frequent_histologicals,
    plot_frequent_topographies,
    plot_frequent_morphologies,
    ncol = 4, nrow = 1),
  
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

