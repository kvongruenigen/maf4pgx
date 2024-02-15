################################################################################
#                         EDA of acquired data                                 #
################################################################################

# Notes ------------------------------------------------------------------------

# Idea:
# Describe the data and analyze what is possible to with the data

# To Dos:
# - Affected genes
# - Gene length normalization
# - Distinguish CNV DEL from SNV DEL

# There are 3 data sets:
# - metadata (individuals + biosamples)
# - cnv variants
# - snv variants

# assign(".Last",  function() {
#                              system("R")}, envir = globalenv())
# quit(save = "no")
################################################################################
#                         Let's start                                          #
################################################################################

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

# Themes for graphs
mytheme <-   theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14))

poster_theme <-   theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 28),
        plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 26),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 26))

# Study population -------------------------------------------------------------

# Load data
metadata <- read_csv("progenetix/progenetix_tcga_individuals_biosamples_combined_data.csv")

# Sort columns in alphabetical order
metadata <- metadata[, order(names(metadata))]

# Convert multiple columns to factors
metadata <- metadata %>%
  mutate_at(vars(stage, substage, histological_diagnosis, icdo_morphology,
                 icdo_topography, sample_origin, tumor_type, project,
                 vital_status, sex, ethnicity),
            as.factor)

# Variables:
dput(names(metadata))
# age_at_diagnosis
# biosample_id
# birth_year
# callset_ids
# collection_moment_years
# days_to_death
# disease
# ethnicity
# followup_time_months
# histological_diagnosis
# icdo_morphology
# icdo_topography
# individual_id
# project
# sample_origin
# sex
# stage
# substage
# tumor_type
# updated
# vital_status

summary(metadata)


## Table 1 ---------------------------------------------------------------------
table1 <- metadata %>%
  select(c(sex, birth_year, age_at_diagnosis, stage,
           ethnicity, vital_status,
           days_to_death, followup_time_months))

table1 <- CreateTableOne(data = table1)
## Why is days_to_death different from survival_time_days?
## survival_time_days is derived from days_to_death -> take days_to_death

ptable1 <- print(table1, printToggle = FALSE, noSpaces = TRUE)
kable(ptable1, format = "latex")


## Table 2 ---------------------------------------------------------------------
table2 <- metadata %>%
  select(c(substage, histological_diagnosis, icdo_morphology,
           icdo_topography, sample_origin, tumor_type, collection_moment_years))
## sample_origin, icdo_topography, icdo_morphology,
## histologicals_diagnosis too long

table2 <- metadata %>%
  select(c(substage, tumor_type, collection_moment_years, project))

CreateTableOne(data = table2)


## Table -----------------------------------------------------------------------
table <- metadata %>%
  select(c(sex, birth_year, age_at_diagnosis, stage,
           ethnicity, vital_status,
           days_to_death, followup_time_months, project))

table1 <- CreateTableOne(data = table)

ptable1 <- print(table1, printToggle = FALSE, noSpaces = TRUE)
kable(ptable1, format = "latex")
# table2 <- CreateTableOne(data = table, strata = c("project"))


## Barcharts -------------------------------------------------------------------



### Stage ----------------------------------------------------------------------

# Stage Bar Chart
ggplot(metadata, aes(stage)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "Stage", x = "", y = "Count") +
  mytheme

# Substage Table
substages <- metadata %>%
  select(c(stage, substage)) %>%
  filter(stage != "Stage Unknown", stage != "Stage X",
         stage != "Stage 0", stage != "Stage I/II NOS")

# Substage Bar Chart --- Don't like the legend
ggplot(substages, aes(stage)) +
  geom_bar(aes(fill = substage), col = "black") +
  labs(title = "Substage", x = "", y = "Count") +
  mytheme

### Sample site ----------------------------------------------------------------
# TBC
ggplot(metadata, aes(sample_origin)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(x = "", y = "Count", title = "Sample sites") +
  mytheme

unique(metadata$sample_origin)

### Tumor type -----------------------------------------------------------------
levels(metadata$tumor_type) <- c("New Primary", "Metastatic", "Primary Blood",
                                  "Primary Tumor", "Recurrent Tumor")
metadata$tumor_type <- factor(
  metadata$tumor_type,
  levels = names(sort(table(metadata$tumor_type), decreasing = TRUE))
)

ggplot(metadata, aes(tumor_type)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "Tumor Type", x = "", y = "Count") +
  scale_y_log10() +
  mytheme

### Project --------------------------------------------------------------------
# # Remove "TCGA" and "project" from project names
# metadata$project <- gsub("TCGA ", "", metadata$project)
# metadata$project <- gsub(" project", "", metadata$project)

# Filter the unique individuals
unique_individuals <- metadata[!duplicated(metadata$individual_id), ]

# Sort the projects by number of individuals
individuals_project <- factor(
  unique_individuals$project,
  levels = names(sort(table(metadata$project), decreasing = FALSE))
)

# Plot the sorted projects
ggplot(unique_individuals, aes(individuals_project)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "TCGA projects", x = "", y = "Individuals") +
  coord_flip() +
  mytheme

table(unique_individuals$project)
unique(metadata$project)

# ct_legend <- list(
#   "ACC" = "Adrenocortical carcinoma",
#   "BLCA" = "Bladder Urothelial Carcinoma",
#   "BRCA" = "Breast invasive carcinoma",
#   "CESC" = "Cervical squamous cell carcinoma and endocervical adenocarcinoma",
#   "CHOL" = "Cholangiocarcinoma",
#   "COAD" = "Colon adenocarcinoma",
#   "DLBC" = "Diffuse large B cell lymphoma",
#   "ESCA" = "Esophageal carcinoma",
#   "GBM" = "Glioblastoma multiforme",
#   "HNSC" = "Head and Neck squamous cell carcinoma",
#   "KICH" = "Kidney Chromophobe",
#   "KIRC" = "Kidney renal clear cell carcinoma",
#   "KIRP" = "Kidney renal papillary cell carcinoma",
#   "LAML" = "Acute Myeloid Leukemia",
#   "LGG" = "Brain Lower Grade Glioma",
#   "LIHC" = "Liver hepatocellular carcinoma",
#   "LUAD" = "Lung adenocarcinoma",
#   "LUSC" = "Lung squamous cell carcinoma",
#   "MESO" = "Mesothelioma",
#   "OV" = "Ovarian serous cystadenocarcinoma",
#   "PAAD" = "Pancreatic adenocarcinoma",
#   "PCPG" = "Pheochromocytoma and Paraganglioma",
#   "PRAD" = "Prostate adenocarcinoma",
#   "READ" = "Rectum adenocarcinoma",
#   "SARC" = "Sarcoma",
#   "SKCM" = "Skin Cutaneous Melanoma",
#   "STAD" = "Stomach adenocarcinoma",
#   "TGCT" = "Testicular Germ Cell Tumors",
#   "THCA" = "Thyroid carcinoma",
#   "THYM" = "Thymoma",
#   "UCEC" = "Uterine Corpus Endometrial Carcinoma",
#   "UCS" = "Uterine Carcinosarcoma",
#   "UVM" = "Uveal Melanoma"
# )

# CNV variants -----------------------------------------------------------------
cnvs <- read_csv("progenetix/progenetix_tcga_cnv_variants_data_with_genes.csv")

cnvs <- cnvs[, order(names(cnvs))]

# Convert multiple columns to factors
cnvs <- cnvs %>%
  mutate_at(vars(chromosome, cnv_state, relative_copy_class, variant_state),
            as.factor)

# Variables:
dput(names(cnvs))
# "affected_genes"
# "biosample_id"
# "callset_id"
# "chromosome"
# "cnv_state"
# "cnv_value"
# "cnvfraction"
# "delfraction"
# "dupfraction"
# "end"
# "individual_id"
# "relative_copy_class"
# "start"
# "updated"
# "variant_id"
# "variant_internal_id"
# "variant_length"
# "variant_state"

summary(cnvs)

# Sequence Alterations ---------------------------------------------------------
snvs <- read_csv("maf_analysis_data.csv")

# Set all columns to lower case
names(snvs) <- tolower(names(snvs))

# Rename column Hugo_Symbol to gene
names(snvs)[names(snvs) == "hugo_symbol"] <- "gene"

# Order columns alphabetically
snvs <- snvs[, order(names(snvs))]

# Convert multiple columns to factors
snvs <- snvs %>%
  mutate_at(vars(biotype, canonical, center, chromosome, clin_sig, gene, hotspot,
                 impact, mirna, one_consequence, variant_classification,
                 variant_type),
            as.factor)

# Variables:
dput(names(snvs))

summary(snvs)

# TO-DOS:
# - How to handle list information in R?
# - aminoacid_changes: Extract list info
# - clinical_interpretations: Extract list info
# - consequence: Extract list info


# Aminoacid changes
unique(snvs$consequence)

# Story telling ----------------------------------------------------------------

# Sorted projects
ggplot(unique_individuals, aes(individuals_project)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "TCGA projects", x = "", y = "Individuals") +
  coord_flip() +
  mytheme

