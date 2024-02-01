################################################################################
#                         EDA of acquired data                                 #
################################################################################

# Notes ------------------------------------------------------------------------

# Idea:
# Describe the data and analyze what is possible to with the data

# To Dos:
# - Affected genes
# - Gene length normalization
# - ClinVar data for CNVs
# - PROBLEM!! INS and DEL are not handled correctly
#   - INS are missing sequence
#   - DEL are missing completly

# Missing in database:
# - Variant Type
# - Variant Classification
# - Affected genes already in MAF data

# Possible from MAF:
# - Variant Type
# - Variant Classification
# - SYMBOL (Affected genes)
# - IMPACT
# - CANONICAL
# - BIOTYPE
# - Consequence
# (- CLIN_SIG)
# (- SIFT)
# (- PolyPhen)


# There are 4 data sets:
# - individuals
# - biosamples
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

# Individuals ------------------------------------------------------------------
individuals <- read_csv("progenetix/progenetix_tcga_individuals_data.csv")

# Variables:
dput(names(individuals))
# updated
# individual_id
# days_to_death
# birthyear
# vital_status
# sex
# external_references
# disease
# stage
# survival_time_days
# ethnicity
# age_at_diagnosis
# followup_time_months

# Converting strings to factors
individuals$vital_status <- as.factor(individuals$vital_status)
individuals$sex <- as.factor(individuals$sex)
individuals$ethnicity <- as.factor(individuals$ethnicity)
## Stage from biosamples is more detailed
individuals$stage <- as.factor(individuals$stage)

summary(individuals)

## Table 1 ---------------------------------------------------------------------
table1 <- individuals %>%
  select(c(sex, birthyear, age_at_diagnosis, stage,
           ethnicity, vital_status,
           days_to_death, followup_time_months))

table1 <- CreateTableOne(data = table1)
## Why is days_to_death different from survival_time_days?
## survival_time_days is derived from days_to_death -> take days_to_death

# Table to latex format
ptable1 <- print(table1, printToggle = FALSE, noSpaces = TRUE)
kable(ptable1, format = "latex")

# Biosamples -------------------------------------------------------------------
biosamples <- read_csv("progenetix/progenetix_tcga_biosamples_data.csv")

# Variables:
dput(names(biosamples))
# "individual_id"
# "substage"
# "callset_ids"
# "biosample_id"
# "histological_diagnosis"
# "icdo_morphology"
# "icdo_topography"
# "tumor_type"
# "sample_origin"
# "collection_moment_years"
# "project"
# "stage"

# Converting strings to factors
biosamples$substage <- as.factor(biosamples$substage)
biosamples$stage <- as.factor(biosamples$stage)
biosamples$histological_diagnosis <- as.factor(biosamples$histological_diagnosis) # nolint: line_length_linter.
biosamples$icdo_morphology <- as.factor(biosamples$icdo_morphology)
biosamples$icdo_topography <- as.factor(biosamples$icdo_topography)
biosamples$sample_origin <- as.factor(biosamples$sample_origin)
biosamples$tumor_type <- as.factor(biosamples$tumor_type)
biosamples$project <- as.factor(biosamples$project)

summary(biosamples)

## Table 2 ---------------------------------------------------------------------
table2 <- biosamples %>%
  select(c(substage, histological_diagnosis, icdo_morphology,
           icdo_topography, sample_origin, tumor_type, collection_moment_years))
## sample_origin, icdo_topography, icdo_morphology,
## histologicals_diagnosis too long

table2 <- biosamples %>%
  select(c(substage, tumor_type, collection_moment_years, project))

CreateTableOne(data = table2)

# Biosamples + Individuals -----------------------------------------------------
individuals <- individuals %>%
  select(-c(stage))
meta_data <- merge(biosamples, individuals, by = "individual_id")

## Table -----------------------------------------------------------------------
table <- meta_data %>%
  select(c(sex, birthyear, age_at_diagnosis, stage,
           ethnicity, vital_status,
           days_to_death, followup_time_months, project))

table1 <- CreateTableOne(data = table)
## Why is days_to_death different from survival_time_days?
## survival_time_days is derived from days_to_death -> take days_to_death

ptable1 <- print(table1, printToggle = FALSE, noSpaces = TRUE)
kable(ptable1, format = "latex")
# table2 <- CreateTableOne(data = table, strata = c("project"))
## Barcharts -------------------------------------------------------------------

### Stage ----------------------------------------------------------------------
ggplot(meta_data, aes(stage)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "Stage", x = "", y = "Count") +
  coord_flip() +
  mytheme

substages <- meta_data %>%
  select(c(stage, substage)) %>%
  filter(stage != "Stage Unknown", stage != "Stage X",
         stage != "Stage 0", stage != "Stage I/II NOS")

ggplot(substages, aes(stage)) +
  geom_bar(aes(fill = substage), col = "black") +
  labs(title = "Substage", x = "", y = "Count") +
  coord_flip() +
  mytheme

### Sample site ----------------------------------------------------------------
# TBC
ggplot(meta_data, aes(sample_origin)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(x = "", y = "Count", title = "Sample sites") +
  mytheme

unique(meta_data$sample_origin)

### Tumor type -----------------------------------------------------------------
levels(meta_data$tumor_type) <- c("New Primary", "Metastatic", "Primary Blood",
                                  "Primary Tumor", "Recurrent Tumor")
meta_data$tumor_type <- factor(
  meta_data$tumor_type,
  levels = names(sort(table(meta_data$tumor_type), decreasing = TRUE))
)

table(meta_data$tumor_type)

ggplot(meta_data, aes(tumor_type)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "Tumor Type", x = "", y = "Samples (log)") +
  scale_y_log10() +
  mytheme

### Project --------------------------------------------------------------------
# Remove "TCGA" and "project" from project names
meta_data$project <- gsub("TCGA ", "", meta_data$project)
meta_data$project <- gsub(" project", "", meta_data$project)

#Filter the unique individuals
unique_individuals <- meta_data[!duplicated(meta_data$individual_id), ]

# Sort the projects by number of individuals
individuals_project <- factor(
  unique_individuals$project,
  levels = names(sort(table(meta_data$project), decreasing = FALSE))
)

# Plot
project_barchart <- ggplot(meta_data, aes(project)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "TCGA projects", x = "", y = "Individuals") +
  coord_flip() +
  mytheme

table(unique_individuals$project)

# CNV variants -----------------------------------------------------------------
cnvs <- read_csv("progenetix/progenetix_tcga_cnv_variants_data_with_genes.csv")

# Variables:
dput(names(cnvs))

sumary(cnvs)

# Sequence Alterations ---------------------------------------------------------
snv_variants <- read_csv("data/matching_maf_data_curated.csv")

# Story telling ----------------------------------------------------------------
ggplot(meta_data, aes(project)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "TCGA projects", x = "", y = "Individuals") +
  coord_flip() +
  mytheme
unique(meta_data$project)
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
