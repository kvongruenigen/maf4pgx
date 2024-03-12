#___________________________________________####################################
# Sequence Alterations ---------------------------------------------------------
#___________________________________________####################################

#"@To-Dos:
# - Extract clinical_information
# - Duplicated values in clin_sig
# - NA in projects



# Variables:
dput(names(snvs))

# Summary
str(snvs)

# Explore numerical variables
numeric_summary_snvs <- snvs %>%
  select_if(is.numeric) %>%
  summary()
print(numeric_summary_snvs)

# Explore categorical variables
categorical_summary_snvs <- snvs %>%
  select_if(is.factor) %>%
  summary()
print(categorical_summary_snvs)

# Split the strings into lists
snvs <- snvs %>%
  mutate(clinvar_effect = strsplit(gsub("[{}']", "", clinvar_effects), ", ")) %>%
  select(-c(clinvar_effects))

snvs <- snvs %>%
  mutate(clinvar_interpretation = strsplit(gsub("[{}']", "", clinvar_interpretations), ", ")) %>%
  select(-c(clinvar_interpretations))


## Visual exploration -----------------------------------------------------------

### Biotype --------------------------------------------------------------------
plot_biotype <-
plot_top_levels(snvs, "biotype", top_n = 30, decreasing = F)+
  scale_y_log10(labels = label_log(digits = 2))+
  labs(title = "Transcript Biotypes",
       y= "Variants (log10)")+
  coord_flip()+
  theme(panel.grid.major.x = element_line(color = "gray", linetype = "solid"))
plot_biotype

### Canonical ------------------------------------------------------------------
# Plot distribution of canonical
plot_canonical <-
ggplot(snvs, aes(x = canonical)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(
    title = "VEP-based Canonical Transcript",
    x = "",
    y = "Variants (log10)") +
  scale_y_log10(labels = label_log(digits = 2)) +
  scale_x_discrete(
    labels = c(
      "YES" = "Yes",
      "NO" = "No")
  ) +
  geom_text(stat = "count", aes(label = ..count..), hjust = 1.5) +
  geom_text(stat = "count", aes(label = scales::percent(..count../sum(..count..))), hjust = -0.2) +
  theme_basic +
  coord_flip()+
  theme(panel.grid.major.x = element_line(color = "gray", linetype = "solid"))
plot_canonical

### Centers --------------------------------------------------------------------
# Plot sorted distribution of centers
snvs$center <- factor(
  snvs$center,
  levels = names(sort(table(snvs$center), decreasing = FALSE))
)
levels(snvs$center)

plot_centers <-
ggplot(snvs, aes(x = center)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "Sequencing Centers", x = "", y = "Variants (log10)") +
  coord_flip() +
  scale_y_log10(labels = label_log(digits = 2)) +
  theme_basic +
  theme(panel.grid.major.x = element_line(color = "gray", linetype = "solid")) +
  geom_text(stat = "count", aes(label = ..count..), hjust = 1.5) +
  geom_text(stat = "count", aes(label = scales::percent(..count../sum(..count..))), hjust = -0.2)
plot_centers

### ClinVar Significance -------------------------------------------------------
# Plot ClinVar significance distribution
levels(snvs$clin_sig)
# Split the entries in the clin_sig column
clin_sig_lists <- strsplit(as.character(snvs$clin_sig), ";")

# Count the occurrences of each term
clin_sig_counts <- table(unlist(clin_sig_lists))

# Convert the counts to a data frame
clin_sig_df <- as.data.frame(clin_sig_counts)
clin_sig_df <- clin_sig_df[order(clin_sig_df$Freq), ]  # Order by frequency

# Convert Var1 to a factor with the correct order
clin_sig_df$Var1 <- factor(clin_sig_df$Var1, levels = clin_sig_df$Var1)

# Create a bar plot
plot_clin_sig <-
ggplot(clin_sig_df, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "lightblue", col="black") +
  labs(title = "Clinical Significance", x = "", y = "Variants (log10)") +
  scale_x_discrete(
    labels = c(
      "benign" = "Benign",
      "likely benign" = "Likely Benign",
      "uncertain significance" = "Uncertain Significance",
      "likely pathogenic" = "Likely Pathogenic",
      "pathogenic" = "Pathogenic",
      "drug response" = "Drug Response",
      "conflicting interpretations of pathogenicity" = "Conflicting Interpretations of Pathogenicity",
      "other" = "Other",
      "risk factor" = "Risk Factor",
      "protective" = "Protective",
      "association" = "Association",
      "not provided" = "Not Provided",
      "benign/likely benign" = "Benign/Likely Benign",
      "affects" = "Affects")
  ) +
  scale_y_log10(labels = label_log(digits = 2)) +
  theme_basic +
  theme(panel.grid.major.x = element_line(color = "gray", linetype = "solid"))+
  coord_flip()
plot_clin_sig

### Consequence Types ----------------------------------------------------------

# Plot the top 3 sample origins - works
plot_consequences <- 
  plot_top_levels(snvs, "one_consequence", top_n = 2) +
  scale_x_discrete(
    labels = c(
      "missense variant" = "Missense",
      "synonymous variant" = "Synonymous"
  ))+
  scale_y_continuous(labels = scales::comma_format(), limits = c(0,2000000)) +
  labs(title = "Most Frequent Consequence Types",
       x = "",
       y = "Variants") +
  theme_basic +
  coord_flip() +
  theme(panel.grid.major.x = element_line(color = "gray", linetype = "solid")) +
  geom_text(stat = "count", aes(label = ..count..), hjust = 1.5) +
  geom_text(stat = "count", aes(label = scales::percent(..count../sum(..count..))), hjust = -0.2)
plot_consequences

# CONSEQUENCES
# Split the entries in the clin_sig column
snvs_filtered <- snvs[!is.na(snvs$consequence), ]

consequence_lists <- strsplit(as.character(snvs_filtered$consequence), ",")

# Count the occurrences of each term
consequence_lists <- lapply(consequence_lists, function(levels_split) str_trim(levels_split))
consequence_counts <- table(unlist(consequence_lists))

# Convert the counts to a data frame
consequence_df <- as.data.frame(consequence_counts)
consequence_df <- consequence_df[order(consequence_df$Freq), ]  # Order by frequency

# Convert Var1 to a factor with the correct order
consequence_df$Var1 <- factor(consequence_df$Var1, levels = consequence_df$Var1)
colnames(consequence_df)[1] <- "consequence"
print(consequence_df)

# Create a bar plot
ggplot(consequence_df, aes(x = consequence, y = Freq)) +
  geom_bar(stat = "identity", fill = "lightblue", col="black") +
  labs(title = "Consequence Types", x = "", y = "Variants (log10)") +
  scale_y_log10(labels = label_log(digits = 2)) +
  theme_basic +
  coord_flip()+
  theme(panel.grid.major.x = element_line(color = "gray", linetype = "solid"))+
  geom_text(stat = "identity", aes(label = scales::percent(Freq/sum(Freq),
                                   accuracy = 0.1)), hjust = 1.5)

### Genes ----------------------------------------------------------------------
# Plot the top 10 mutated genes
plot_top_mutated_genes <- 
  plot_top_levels(snvs, "gene", top_n = 20, include_others = F) +
  labs(title = "Top 20 Mutated Genes (Seuqence Variations)",
       y = "Variants per Gene") +
  scale_y_continuous(labels = scales::comma_format()) +
  coord_flip()+
  theme_basic +
  theme(panel.grid.major.x = element_line(color = "gray", linetype = "solid"))
plot_top_mutated_genes

# Show number and percentage for top genes
top_genes <- snvs %>%
  count(gene) %>%
  arrange(desc(n)) %>%
  top_n(20) %>%
  mutate(percentage = n / sum(n))
top_genes

# Show number and percentage for top genes
top_genes_project <- snvs %>%
  group_by(project) %>%
  count(gene) %>%
  arrange(desc(n)) %>%
  top_n(20) %>%
  mutate(percentage = n / sum(n))
top_genes_project

# Exclude projects with 100 mutations per gene
top_genes_project_filtered <- top_genes_project %>%
  filter(n > 200)

# Show top genes per project
ggplot(top_genes_project_filtered, aes(x = reorder(gene, n), y = n, fill = project)) +
  geom_bar(stat = "identity", col = 'black') +
  labs(title = "Top 20 Mutated Genes per Project",
       x = "Gene",
       y = "Variants") +
  scale_y_continuous(labels = scales::comma_format()) +
  coord_flip() +
  theme_classic() +
  theme(panel.grid.major.x = element_line(color = "gray", linetype = "solid")) +
  theme(legend.position = "bottom")


### Hotspots -------------------------------------------------------------------
# Plot distribution of hotspots
plot_hotspots <-
ggplot(snvs, aes(x = hotspot)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(
    title = "Mutational Hotspot",
    x = "Hotspot",
    y = "Variants (log10)") +
  scale_y_log10(labels = label_log(digits = 2)) +
  scale_x_discrete(
    labels = c(
      "Y" = "Hotspot",
      "N" = "Not Hotspot")
  ) +
  theme_basic +
  coord_flip()+
  theme(panel.grid.major.x = element_line(color = "gray", linetype = "solid"))+
  geom_text(stat = "count", aes(label = ..count..), hjust = 1.5) +
  geom_text(stat = "count", aes(label = scales::percent(..count../sum(..count..))), hjust = -0.2)

plot_hotspots


### Impact of Variants --------------------------------------------------------

# Sort
snvs$impact_sorted <- factor(
  snvs$impact,
  levels = names(sort(table(snvs$impact), decreasing = FALSE))
)
table(snvs$impact)
# Plot
plot_impact <- 
  ggplot(snvs, aes(x = factor(impact_sorted))) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(
    title = "Impact Modifier for Consequence",
    x = "",
    y = "Variants (log10)") +
  scale_x_discrete(
    labels = c(
      "LOW" = "Low",
      "MODERATE" = "Moderate",
      "HIGH" = "High",
      "MODIFIER" = "Modifier")
    ) +
  scale_y_log10(labels = label_log(digits = 2), limits = c(1,10^7)) +
  theme_basic +
  coord_flip()+
  theme(panel.grid.major.x = element_line(color = "gray", linetype = "solid"))+
  geom_text(stat = "count", aes(label = ..count..), hjust = 1.5) +
  geom_text(stat = "count", aes(label = scales::percent(..count../sum(..count..))), hjust = -0.2)

plot_impact


### miRNA ---------------------------------------------------------------------
# Split the entries in the miRNA column
miRNA_lists <- strsplit(as.character(snvs$mirna), ";")

sum(!is.na(snvs$mirna))

# Count the occurrences of each term
miRNA_counts <- table(unlist(miRNA_lists))

# Convert the counts to a data frame
miRNA_df <- as.data.frame(miRNA_counts)
miRNA_df <- miRNA_df[order(miRNA_df$Freq), ]  # Order by frequency

# Convert Var1 to a factor with the correct order
miRNA_df$Var1 <- factor(miRNA_df$Var1, levels = miRNA_df$Var1)

# Create a bar plot
plot_mirna <-
ggplot(miRNA_df, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "lightblue", col="black") +
  labs(title = "miRNA", x = "", y = "Variants (log10)") +
  scale_y_log10(labels = label_log(digits = 2)) +
  theme_basic +
  coord_flip()+
  # Add percentage and absolute values
  geom_text(aes(label = Freq), hjust = 1.5) +
  geom_text(aes(label = scales::percent(Freq/sum(Freq))), hjust = -0.2)+
  theme(panel.grid.major.x = element_line(color = "gray", linetype = "solid"))
plot_mirna

### Project: Mutation per project ----------------------------------------------
# Plot the distribution of mutations per project
# Need to adjust by number of samples
snvs$project <- factor(
  snvs$project,
  levels = names(sort(table(snvs$project), decreasing = FALSE))
)
plot_mutations_per_project <- 
  ggplot(snvs, aes(x = project)) +
  geom_bar(fill = "lightblue", col = "black", na.rm = T) +
  labs(
    title = "Mutations per Project",
    x = "Project",
    y = "Variants (log10)") +
  scale_y_log10(labels = label_log(digits = 2)) +
  theme_basic +
  theme(panel.grid.major.x = element_line(color = "gray", linetype = "solid"))+
  coord_flip()
plot_mutations_per_project


# Calculate the number of variants per project
variants_per_project <- snvs %>%
  group_by(project, biosample_id) %>%
  summarise(variant_count = n())
avg_variants <- variants_per_project%>%
  group_by(project) %>%
  summarise(avg_variant_count = mean(variant_count))

# Join the calculated averages back to your main dataset
variants_per_project <- left_join(variants_per_project, avg_variants, by = "project")

# Create a new column for the normalized variant count
variants_per_project$normalized_variant_count <- variants_per_project$variant_count / variants_per_project$avg_variant_count

# Order the levels of the project factor by the average variant count
variants_per_project$project <- factor(variants_per_project$project, levels = levels(factor(variants_per_project$project))[order(variants_per_project$avg_variant_count)])

# Create the plot
plot_mutations_per_project <- 
  ggplot(variants_per_project, aes(x = project, y = normalized_variant_count)) +
  geom_bar(stat = "identity", fill = "lightblue", col = "black", na.rm = TRUE) +
  labs(
    title = "Normalized Mutations per Project",
    x = "Project",
    y = "Normalized Variants (log10)"
  ) +
  scale_y_log10(labels = scales::label_log(base = 10)) +
  theme_minimal() +  # Change the theme if needed
  theme(panel.grid.major.x = element_line(color = "gray", linetype = "solid")) +
  coord_flip()

# Display the plot
print(plot_mutations_per_project)


### PolyPhen -------------------------------------------------------------------
# Sort
snvs$polyphen_sorted <- factor(
  snvs$polyphen,
  levels = names(sort(table(snvs$polyphen), decreasing = FALSE))
)
# Plot distribution of PolyPhen
plot_polyphen <-
ggplot(snvs, aes(x = polyphen)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(
    title = "PolyPhen",
    x = "",
    y = "Variants (log10)") +
  scale_y_log10(labels = label_log(digits = 2)) +
  scale_x_discrete(
    labels = c(
      "benign" = "Benign",
      "possibly_damaging" = "Possibly Damaging",
      "probably_damaging" = "Probably Damaging",
      "unknown" = "Unknown")
  ) +
  theme_basic +
  theme(panel.grid.major.x = element_line(color = "gray", linetype = "solid"))+
  coord_flip()+
  geom_text(stat = "count", aes(label = ..count..), hjust = 1.5) +
  geom_text(stat = "count", aes(label = scales::percent(..count../sum(..count..))), hjust = -0.2)
plot_polyphen

### SIFT -----------------------------------------------------------------------
# Sort
snvs$sift_sorted <- factor(
  snvs$sift,
  levels = names(sort(table(snvs$sift), decreasing = FALSE))
)

snvs_filtered <- snvs[!is.na(snvs$sift), ]



# Plot distribution of SIFT
plot_sift <-
ggplot(snvs_filtered, aes(x = sift)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(
    title = "SIFT",
    x = "",
    y = "Variants (log10)") +
  scale_y_log10(labels = label_log(digits = 2)) +
  scale_x_discrete(
    labels = c(
      "deleterious" = "Deleterious",
      "tolerated" = "Tolerated",
      "tolerated_low_confidence" = "Tolerated (low confidence)",
      "deleterious_low_confidence" = "Deleterious (low confidence)",
      "unknown" = "Unknown")
  ) +
  theme_basic +
  theme(panel.grid.major.x = element_line(color = "gray", linetype = "solid"))+
  coord_flip()+
  geom_text(stat = "count", aes(label = ..count..), hjust = 1.5) +
  geom_text(stat = "count", aes(label = scales::percent(..count../sum(..count..))), hjust = -0.2)
plot_sift

# Together
ggarrange(
  plot_polyphen,
  plot_sift,
  ncol = 2
)

### Interpretations table ------------------------------------------------------
# Count the occurrences of each term for each column
# clin_sig
# sift
# polyphen
interpretations_table <- snvs %>%
  select(clin_sig, sift, polyphen) %>%
  filter(!is.na(clin_sig) & !is.na(sift) & !is.na(polyphen)) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
  mutate(value = strsplit(as.character(value), ";")) %>%
  unnest(value) %>%
  group_by(variable, value) %>%
  summarise(n = n()) %>%
  arrange(variable, desc(n))

### Variation Classification ---------------------------------------------------
# Sort by frequency
snvs$variant_classification_sorted <- factor(
  snvs$variant_classification,
  levels = names(sort(table(snvs$variant_classification), decreasing = FALSE))
)

# Make labels pretty
levels(snvs$variant_classification_sorted) <- gsub("_", " ", levels(snvs$variant_classification_sorted))


# Plot
plot_variant_classification <-
  ggplot(snvs, aes(x = variant_classification_sorted)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "Variant Classification", x = "", y = "Variants (log10)") +
  scale_y_log10(labels = label_log(digits = 2)) +
  coord_flip()+
  theme_basic+
  theme(panel.grid.major.x = element_line(color = "gray", linetype = "solid"))
plot_variant_classification

### Variation Types ------------------------------------------------------------
# Sort by frequency
snvs$variant_type_sorted <- factor(
  snvs$variant_type,
  levels = names(sort(table(snvs$variant_type), decreasing = FALSE))
)

# Plot distribution of variant types
plot_variant_types <- 
  ggplot(snvs, aes(x = variant_type_sorted)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "Variant Types", x = "", y = "Variants (log10)")+
  scale_y_log10(labels = label_log(digits = 2)) +
  theme_basic +
  coord_flip() +
  theme(panel.grid.major.x = element_line(color = "gray", linetype = "solid"))+
  geom_text(stat = "count", aes(label = ..count..), hjust = 1.5) +
  geom_text(stat = "count", aes(label = scales::percent(..count../sum(..count..))), hjust = -0.2)
plot_variant_types

### Variations per sample and project ------------------------------------------

variations_per_biosample_project <- snvs %>%
  group_by(biosample_id, project) %>%
  filter(!is.na(project)) %>% # Remove NA values
  summarise(num_variants = n(), cnvfraction) %>%
  distinct(biosample_id, project, .keep_all = TRUE)

plot_variations_per_sample <- 
  ggplot(variations_per_biosample_project,
         aes(x = reorder(project, -num_variants, median), y = num_variants)) +
  geom_violin(fill = "lightblue", width = 1.23, trim = F) +
  geom_boxplot(width = 0.07, outlier.size = 0, col = "black") +
  scale_y_log10(labels = label_log(digits = 2)) +
  labs(title = "TCGA Program (Tumor Type)",
       x = "",
       y = "Variations per sample (log10)") +
  theme_basic


# Selection of different distributions
vpbp <- variations_per_biosample_project %>%
  filter(project == "BLCA" | project == "STAD")

plot_variations_per_sample_different_distr <-
  ggplot(vpbp, aes(x = reorder(project, -num_variants, median), y = num_variants)) +
  geom_violin(
    fill = "lightblue",
    #fill = "#31a0ff",
    alpha = 0.7) +
  geom_boxplot(width=0.07, outlier.size = 0, col = "navy") +
  scale_y_log10(labels = label_log(digits = 2))+
  labs(title = "Example of Different Distributions",
       x = "Cancer type",
       y = "Number of Variants (log10)") +
  theme_basic
  
plot_variations_per_sample_different_distr


### Variants vs CNV fraction ---------------------------------------------------
## Per project
plot_snv_cnv_fraction <-
  ggplot(variations_per_biosample_project, aes(log(num_variants), cnvfraction)) +
  geom_point(aes(col = project), show.legend = FALSE) +
  stat_smooth(method = "lm", formula = y ~ x, geom = "smooth") +
  facet_wrap(~project)+
  labs(
    title = "CNV Fraction vs Sequence Variations",
    subtitle = "Per Project",
    x = "Variants (log10)",
    y = "CNV Fraction")+
  theme_basic +
  theme(panel.grid.major.x = element_line(color = "gray", linetype = "dashed"),
        panel.grid.major.y = element_line(color = "gray", linetype = "dashed"))

plot_snv_cnv_fraction


# All plots together -----------------------------------------------------------
plots_snvs_complete <-
ggarrange(
  plot_biotype,
  plot_canonical,
  plot_centers,
  plot_clin_sig,
  plot_consequences,
  plot_top_mutated_genes,
  plot_hotspots,
  plot_impact,
  plot_mirna,
  plot_mutations_per_project,
  plot_polyphen,
  plot_sift,
  plot_variant_classification,
  plot_variant_types,
  plot_variations_per_sample,
  plot_variations_per_sample_different_distr,
  plot_snv_cnv_fraction,
  ncol = 4, nrow = 5
)
plots_snvs_complete


