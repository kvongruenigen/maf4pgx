#___________________________________________####################################
# CNV variants -----------------------------------------------------------------
#___________________________________________####################################

# Variables:
dput(names(cnvs))
# affected_genes          
# biosample_id            #
# callset_id              #
# chromosome              #
# cnv_value               #
# cnvfraction             #
# delfraction             #
# dupfraction             #
# end                     #
# individual_id           #
# start                   #
# updated                 #
# variant_id              #
# variant_internal_id     #
# variant_length          #
# variant_state           #


# Number of CNVs
dim(cnvs) # 1833533 CNVS


# Summary
summary(cnvs)


# Explore numerical variables
numeric_summary <- cnvs %>%
  select_if(is.numeric) %>%
  select(-c(end, start)) %>% # Remove start and end
  summary()
print(numeric_summary)


# Explore categorical variables
# categorical_summary <- cnvs %>%
#   select_if(is.factor) %>%
#   summary()
# print(categorical_summary)

# Actually, variant_state is the only valuable categorical variable
summary(cnvs$variant_state)


## Visual exploration ----------------------------------------------------------

### CNV state ------------------------------------------------------------------
plot_cnv_state <-
  ggplot(cnvs, aes(variant_state))+
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "Variant State", x = "", y = "Variations") +
  scale_x_discrete(labels = c(
    "copy number loss" = "Copy Number Loss",
    "copy number gain" = "Copy Number Gain"
  )) +
  scale_y_continuous(labels = scales::comma_format()) +
  theme_basic +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "solid"))
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5)+
  geom_text(stat = "count", aes(label = scales::percent(..count../sum(..count..))), vjust = 1.5)


### CNV Length -----------------------------------------------------------------

plot_cnv_length <-
  ggplot(cnvs, aes(variant_length)) +
  geom_histogram(fill = "lightblue", col = "black") +
  labs(title = "CNV Length", x = "Nucleotides", y = "Variants (log10)") +
  scale_x_continuous(labels = scales::comma_format())+
  # Make y log scale with 10^x notation
  scale_y_log10(labels = label_log(digits = 2)) +
  theme_basic +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "solid"))
# 0-Inflation -> distribution = Poisson / quasi-poisson?


### CNV Values -------------------------------------------------
plot_cnv_values <-
  ggplot(cnvs, aes(x = cnv_value)) +
  geom_histogram(binwidth = 0.1, fill = "lightblue", col = "black") +
  scale_y_log10(labels = label_log(digits = 2)) +
  labs(title = "CNV Values", x = "", y = "Variants (log10)") +
  theme_basic +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "solid"))


### Fraction Information --------------------------------------------------------
plot_cnv_fractions <-
  ggarrange(
    ggplot(cnvs, aes(x = cnvfraction)) +
      geom_histogram(binwidth = 0.05, fill = "lightblue", col = "black") +
      labs(title = "Distribution of CNV Fraction", x = "Overall Fraction", y = "Frequency") +
      theme_basic,
    ggarrange(
      ggplot(cnvs, aes(x = delfraction)) +
        geom_histogram(binwidth = 0.05, fill = "lightblue", col = "black") +
        labs(title = "", x = "Deletion Fraction", y = "Frequency") +
        theme_basic +
        theme(panel.grid.major.y = element_line(color = "gray", linetype = "solid")),
      
      ggplot(cnvs, aes(x = dupfraction)) +
        geom_histogram(binwidth = 0.05, fill = "lightblue", col = "black") +
        labs(title = "", x = "Duplication Fraction", y = "Frequency") +
        theme_basic +
        theme(panel.grid.major.y = element_line(color = "gray", linetype = "solid")),
      ncol = 2, align = "h"),
    nrow = 2, heights = c(1, 1))
# Normally distributed

### CNV Plots ------------------------------------------------------------------
plots_cnv_complete <-
ggarrange(
  plot_cnv_state,
  ggarrange(
    plot_cnv_fractions,
    ncol=1, nrow=1),
  plot_cnv_length,
  plot_cnv_values,
  ncol = 1, nrow = 4
)
plots_cnv_complete
### Genes ----------------------------------------------------------------------
class(cnvs$affected_genes)

# Remove brackets, " " and "'" from affected_genes
cnvs$affected_genes <- gsub("\\[|\\]|\\'|\\ ", "", cnvs$affected_genes)

# Store affected genes as list in 'affected_genes' column
cnvs$affected_genes <- strsplit(cnvs$affected_genes, ",")


# CNVs that do not affect a single gene
dim(cnvs %>%
      filter(affected_genes == ""))[1]
# 295,933 CNVs do not affect a single gene

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Use separate_rows to split the genes in the 'affected_genes' column
# Then use count to count the occurrences of each gene
# gene_counts <- cnvs %>%
#   separate_rows(affected_genes, sep = ",") %>%
#   count(affected_genes, name = "cnv_count")


# View the resulting data frame with sorted gene counts
# gene_counts <- gene_counts %>%
#   arrange(desc(cnv_count))
# 
# print(gene_counts)
