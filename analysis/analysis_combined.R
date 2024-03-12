#___________________________________________####################################
# Combined ---------------------------------------------------------------------
#___________________________________________####################################

# Combine the data sets --------------------------------------------------------

# Mappings
# biosample/individual id
# location: chromosome + start + end
combine_cnvs <- cnvs %>%
  select(-c(callset_id, individual_id, updated,
            variant_internal_id, variant_id, project))  %>%
  rename(cnv_start = start, cnv_end = end)

combine_snvs <- snvs %>%
  select(-c(project, cnvfraction)) %>%
  rename(snv_start = start, snv_end = end)



# Step 1: Merge metadata with snvs based on biosample_id or individual_id
merged_data <- combine_snvs %>%
  left_join(metadata, by = c("biosample_id", "individual_id"))
str(merged_data)


# Step 2: Merge merged_data with cnvs based on biosample_id and cnv_chromosome
final_merged_data <- merged_data %>%
  # Join over biosample_id and chromosome
  left_join(combine_cnvs, by = c("biosample_id", "chromosome"),
            relationship = "many-to-many") %>%
  # Filter for overlapping regions
  filter(snv_start <= cnv_end & snv_end >= cnv_start)

print(final_merged_data)
str(final_merged_data)


# Genes ------------------------------------------------------------------------
## Gene information -------------------------------------------------------------
str(genes)
combine_genes <- genes %>%
  select(c(gene, gene_end, gene_start, chromosome, gene_length, tier))

complete_data <- final_merged_data %>% 
  # Exlcude all columns containing "sorted"
  select(-contains("sorted")) %>%
  # Join over gene and chromosome
  left_join(combine_genes, by = c("gene", "chromosome"),
            relationship = "many-to-many")
  


# Sort columns alphabetically
complete_data <- complete_data %>%
  select(sort(colnames(complete_data)))

levels(complete_data$variant_state) <- c("DUP", "DEL")

str(complete_data)


## Overall ---------------------------------------------------------------------
mutations_per_gene <- complete_data %>%
  group_by(gene) %>%
  summarize(variants = n()) %>%
  arrange(desc(variants))

mutations_per_gene <- mutations_per_gene %>%
  left_join(combine_genes, by = "gene") %>%
  mutate(variants_per_kb = variants / gene_length)

# Print the summary table
print(mutations_per_gene)

# print only genes with tier
print(mutations_per_gene %>%
        filter(!is.na(tier)),
      n=20)

# Normalized by gene length
print(mutations_per_gene %>%
        arrange(desc(variants_per_kb)) %>%
        filter(!is.na(tier)),
      n=20)

## Per project -----------------------------------------------------------------
mutations_per_project <- complete_data %>%
  group_by(project, gene) %>%
  summarize(variants = n()) %>%
  left_join(combine_genes, by = "gene") %>%
  mutate(variants_per_kb = variants / gene_length) %>%
  arrange(project, desc(variants))

# Print the summary table
print(mutations_per_project, n=20)


# Keep only genes with tier
mutations_per_project_filtered <- mutations_per_project %>%
  filter(!is.na(tier))

print(mutations_per_project_filtered, n=20)

# Keep only top 20 genes per project
top_mutations_per_project <- mutations_per_project_filtered %>%
  group_by(project) %>%
  top_n(5, variants) %>%
  arrange(project, desc(variants))

print(top_mutations_per_project, n=60)
length(unique(top_mutations_per_project$gene))

plot_top_levels(top_mutations_per_project, "gene", top_n = 20,
                include_others = F)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

################################################################################

tsgs <- complete_data %>%
  filter("TSG" %in% cancer_role)

ogs <- complete_data %>%
  filter("oncogene" %in% cancer_role)

double_role <- complete_data %>%
  filter("TSG" %in% cancer_role & "oncogene" %in% cancer_role)
unique(double_role$gene)



str(complete_data)



# check for duplications in location with same biosample id
complete_data %>%
  group_by(biosample_id, chromosome, snv_start, snv_end) %>%
  summarize(n = n()) %>%
  filter(n > 1)

snvs %>%
  group_by(biosample_id, chromosome, start, end) %>%
  summarize(n = n()) %>%
  filter(n > 1)

## Overlapping genes -----------------------------------------------------------
# Make a list of affected genes per sample for snvs
snv_genes <- snvs %>%
  group_by(biosample_id) %>%
  summarise(snv_genes = list(unique(gene)))

# Make a list of affected genes per sample for cnvs
cnv_genes <- cnvs %>%
  # Unnest the affected_genes list
  unnest(affected_genes) %>%
  # Remove any empty strings or NAs
  filter(affected_genes != "" & !is.na(affected_genes)) %>%
  # Group by biosample_id
  group_by(biosample_id) %>%
  # Summarize the unique genes into a list
  summarise(cnv_genes = list(unique(affected_genes)))

# Load projects
projects <- metadata %>%
  select(biosample_id, project)

# Merge the two data frames based on biosample_id
merged_genes <- snv_genes %>%
  left_join(cnv_genes, by = "biosample_id") %>%
  
  # Make a new column for genes that are in snv_genes and cnv_genes
  mutate(overlap = map2(snv_genes, cnv_genes, intersect)) %>%
  
  # Count overlaps
  mutate(overlap_count = map_int(overlap, length)) %>%
  
  # Sort by overlap_count
  arrange(desc(overlap_count)) %>%
  
  # Add project id
  left_join(projects, by = "biosample_id")

rm(projects)

# # Check if unique list of genes
# length(as.list(merged_genes$snv_genes[[1]])) == length(as.list(merged_genes$snv_genes[[1]]))

# Print the resulting data frame
print(merged_genes)

# Plot distribution of overlaps per project
ggplot(merged_genes, aes(x = overlap_count)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(~project, scales = "free") +
  labs(title = "Distribution of Two-Hit Genes",
       subtitle = "Grouped by Project",
       x = "Number of Overlapping Genes",
       y = "Samples") +
  theme_basic

# Multi-Density plot with project as color
plot_density_of_overlaps <-
ggplot(merged_genes, aes(x = overlap_count, color = project)) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of Two-Hit Genes",
       subtitle = "Grouped by Project",
       x = "Number of Overlapping Genes",
       y = "Samples (density)") +
  scale_x_log10()+
  theme_basic +
  theme(#legend.position = "bottom",
        panel.grid.major.y = element_line(color = "gray", linetype = "dashed"))

################################################################################

# Make a list with all the genes per biosample
all_genes <- merged_genes %>%
  select(biosample_id, snv_genes, cnv_genes, overlap) %>%
  pivot_longer(cols = c(snv_genes, cnv_genes, overlap),
               names_to = "mutation_type",
               values_to = "genes") %>%
  unnest(genes) %>%
  select(biosample_id, genes) %>%
  distinct()  # Remove duplicate rows
cnvgene <- cnvs %>%
  select(biosample_id, affected_genes, project) %>%
  filter(affected_genes != "") %>%
  group_by(biosample_id, project) %>%
  summarise(cnv_genes = list(unique(affected_genes))) %>%
  unnest(cnv_genes)
cnvgene <- cnvs %>%
  select(biosample_id, affected_genes, project) %>%
  filter(affected_genes != "") %>%
  group_by(biosample_id, project) %>%
  summarise(cnv_genes = list(affected_genes)) %>%
  unnest(cnv_genes)

# Unnest the lists of genes
merged_genes_long <- merged_genes %>%
  pivot_longer(cols = c(snv_genes, cnv_genes, overlap),
               names_to = "mutation_type",
               values_to = "genes") %>%
  unnest(genes) %>%
  select(biosample_id, genes, project) %>%
  distinct()  # Remove duplicate rows

# Count occurrences of each gene per project
gene_counts <- merged_genes_long %>%
  group_by(project, genes) %>%
  summarise(mutations = n()) %>%
  arrange(project, desc(mutations))

# Get the top mutated genes per project
top_genes_per_project <- gene_counts %>%
  group_by(project) %>%
  top_n(5, mutations) %>%
  arrange(project, desc(mutations))

# View the result
print(top_genes_per_project)

## Number of variations per gene -----------------------------------------------
### SNVs
snv_gene_counts <- snvs %>%
  count(gene, name = "snv_count")%>%
  arrange(desc(snv_count))

### CNVs
cnv_gene_counts <- cnvs %>%
  #filter(affected_genes != "") %>%
  count(affected_genes, name = "cnv_count") %>%
  arrange(desc(cnv_count))

# Count how many times a gene is mentioned per project
snv_gene_counts_per_project <- snvs %>%
  group_by(project, gene) %>%
  summarise(mutations = n()) %>%
  arrange(project, desc(mutations))


# Count how many times each gene is mutated per sample
cnv_gene_counts_per_project <- cnvgene %>%
  group_by(project, biosample_id, cnv_genes) %>%
  summarise(mutations = n()) %>%
  arrange(project, desc(mutations))



# Double Hits ------------------------------------------------------------------
# Find mutation that overlap within the same sample

str(cnvs)
# biosample_id, chromosome, start, end, 
# cnv_value, cnvfraction, delfraction, dupfraction, variant_state, variant_length

# Map cnvs to snvs
hitting <- snvs %>%
  inner_join(cnvs, by = "biosample_id", relationship = "many-to-many") %>%
  filter(start >= cnv_start & end <= cnv_end)

double_hits <- snvs %>%
  inner_join(cnvs, by = "biosample_id") %>%
  filter(gene.x == gene.y) %>%
  select(biosample_id, gene.x, chromosome.x, start.x, end.x, variant_id.x, variant_id.y)

# View the result
