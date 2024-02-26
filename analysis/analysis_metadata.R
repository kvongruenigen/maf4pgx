#___________________________________________####################################
# Study population -------------------------------------------------------------
#___________________________________________####################################

# Variables:
dput(names(metadata))
# age_at_diagnosis          #
# biosample_id              #
# birth_year                #
# callset_ids               #
# collection_moment_years   
# days_to_death             
# disease                   
# ethnicity                 #
# followup_time_months      #
# histological_diagnosis
# icdo_morphology
# icdo_topography
# individual_id             #
# project                   #
# sample_origin
# sex                       #
# stage                     #
# substage                  #
# tumor_type                #
# updated                   #
# vital_status              #

summary(metadata)


# days_to_death != followup_time_months !!! ------------------------------------
ggplot(metadata, aes(x = followup_time_months, y = days_to_death/30))+
  geom_point()+
  geom_smooth(method = 'lm')+
  theme_minimal()


# Explore the data -------------------------------------------------------------
# Explore numerical variables
numeric_summary <- metadata %>%
  select_if(is.numeric) %>%
  summary()
print(numeric_summary)


# Explore categorical variables
categorical_summary <- metadata %>%
  select_if(is.factor) %>%
  summary()
print(categorical_summary)

# Individuals ----------------------------------------------------------------
## Tables ----------------------------------------------------------------------
table1 <- metadata %>%
  select(c(sex, birth_year, age_at_diagnosis, stage,
           ethnicity, vital_status,
           days_to_death, followup_time_months))

table1 <- CreateTableOne(data = table1)
## Why is days_to_death different from survival_time_days?
## survival_time_days is derived from days_to_death -> take days_to_death

ptable1 <- print(table1, printToggle = FALSE, noSpaces = TRUE)
kable(ptable1, format = "latex")

### Ethnicity
table(metadata$ethnicity)

## Visual exploration ----------------------------------------------------------

### Sex ------------------------------------------------------------------------
table(metadata$sex)

plot_sex <- 
  ggplot(metadata, aes(sex))+
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "Genotypic Sex", x = "", y = "Individuals") +
  scale_x_discrete(labels = c(
    "female genotypic sex" = "Female",
    "male genotypic sex" = "Male"
  )) +
  theme_basic +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "solid")) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5)+
  geom_text(stat = "count", aes(label = scales::percent(round(..count../sum(..count..), 2))), vjust = 1.5)
plot_sex

# Combined with vital status
ggplot(metadata, aes(sex, fill= vital_status))+
  geom_bar(aes(fill = vital_status), col = "black") +
  labs(title = "Sex", x = "", y = "Individuals") +
  theme_basic +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5)+
  geom_text(stat = "count", aes(label = scales::percent(round(..count../sum(..count..), 2))), vjust = 1.5)

### Stage !! -------------------------------------------------------------------

# Stage Bar Chart
plot_stage <-
  ggplot(metadata, aes(stage)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "Stage", x = "", y = "Count") +
  theme_basic +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5)+
  geom_text(stat = "count", aes(label = scales::percent(round(..count../sum(..count..), 2))), vjust = 1.5)


# Substage Bar Chart --- Don't like the legend
ggplot(metadata, aes(stage)) +
  geom_bar(aes(fill = substage), col = "black", na.rm = TRUE) +
  labs(title = "Substage", x = "", y = "Count") +
  theme_basic

### Age distribution -----------------------------------------------------------

plot_age_simple <-
  ggplot(metadata, aes(x = age_at_diagnosis)) +
  geom_histogram(binwidth = 5, fill = "lightblue", col = "black") +
  labs(title = "Age Distribution at Diagnosis", x = "Age", y = "Count") +
  geom_vline(aes(xintercept = mean(age_at_diagnosis, na.rm = TRUE)), size = 1,
             lty = 2, col = 'red', alpha = .8, show.legend = FALSE)+
  theme_basic
# Slight left skew

# Plot age distribution per sex
plot_age_by_sex <-
  ggarrange(
    ggplot(metadata, aes(age_at_diagnosis)) +
      geom_histogram(aes(y =after_stat(density)), fill = 'lightblue',
                     col = 'black',binwidth = 5)+
      geom_vline(aes(xintercept = mean(age_at_diagnosis, na.rm = TRUE)),
                 linewidth = 1, lty = 2, col = 'red', alpha = .8,
                 show.legend = FALSE) +
      labs(title = 'Age distribution',
           subtitle = 'Mean overlaid') +
      xlab('Age [years]') +
      theme_basic,
    
    ggplot(metadata, aes(age_at_diagnosis))+
      geom_histogram(aes(y = ..density..), fill = 'lightblue', col = 'black', binwidth = 5) +
      facet_wrap(~sex) +
      labs(title = '', subtitle = 'Grouped by sex') + 
      xlab('Age [years]') +
      theme_basic,
    ncol = 1, nrow = 2)

plot_age_per_project <-
ggplot(metadata, aes(age_at_diagnosis, fill = sex)) +
  geom_histogram(aes(y =after_stat(density)),
                 binwidth = 5)+
  geom_vline(aes(xintercept = mean(age_at_diagnosis, na.rm = TRUE)),
             linewidth = 0.5, lty = 2, col = 'grey30', alpha = .8,
             show.legend = FALSE) +
  labs(title = 'Age distribution',
       subtitle = 'Overall mean overlaid') +
  facet_wrap(~project) +
  xlab('Age [years]') +
  ylab('Density') +
  theme_basic +
  scale_fill_manual(name = "Genotypic Sex",
                    values = c(my_color_palette[4], my_color_palette[5]), # Change color (F,M)
                    labels = c("female genotypic sex" = "Female",
                               "male genotypic sex" = "Male")) +  # Change legend title and colors
  scale_y_continuous(labels = scales::number_format(scale = 1, accuracy = 0.01))
  

### Ethnicity ------------------------------------------------------------------
# Reorder the levels of race based on frequency
metadata$ethnicity_sorted <- factor(
  metadata$ethnicity,
  levels = names(sort(table(metadata$ethnicity),
                      decreasing = TRUE)
  )
)

levels(metadata$ethnicity_sorted) <- c('white', 'not reported', 'black', 'asian',
                                       'hispanic', 'pacific islander', 'native american')

# Plot ethnicity distribution
plot_ethnicities <-
  ggplot(metadata, aes(x = ethnicity_sorted)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "Ethnicities", x = "", y = "Individuals") +
  scale_x_discrete(labels = c(
    "white" = "White",
    "not reported" = "Not reported",
    "black" = "Black",
    "asian" = "Asian",
    "hispanic" = "Hispanic",
    "pacific islander" = "Pacific Islander",
    "native american" = "Native American"
  )) +
  scale_y_continuous(labels = scales::comma_format()) +
  theme_basic +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "solid")) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5)+
  geom_text(stat = "count", aes(label = scales::percent(round(..count../sum(..count..), 2))), vjust = 1.5)
plot_ethnicities

# Colorful with legend
ggplot(metadata, aes(ethnicity_sorted)) +
  geom_bar(aes(fill = ethnicity_sorted), fill = "lightblue", col = "black") +
  labs(title = 'Ethnicity', fill = '') +
  xlab('Ethnicity') +
  ylab('Amount') +
  theme_basic +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())



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
plot_projects <-
  ggplot(unique_individuals, aes(individuals_project)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "TCGA Projects (Tumor Types)", x = "", y = "Individuals") +
  coord_flip() +
  theme_basic +  
  # geom_text(stat = "count", aes(label = ..count..), hjust = 1.5) +
  theme(panel.grid.major.x = element_line(color = "gray", linetype = "solid"))
plot_projects

table(unique_individuals$project)
unique(metadata$project)

### Follow-up ------------------------------------------------------------------
# Follow-up time
plot_followup_time <-
  ggplot(metadata, aes(followup_time_months))+
  geom_histogram(fill = 'lightblue', col = 'black')+
  geom_vline(aes(xintercept = mean(followup_time_months, na.rm = TRUE)), size = 1,
             lty = 2, col = 'red', alpha = .8, show.legend = FALSE)+
  labs(title = 'Followup time', subtitle = 'Mean overlaid')+
  xlab('Months')+
  ylab('Individuals')+
  scale_y_continuous(labels = scales::comma_format())+
  theme_basic +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "solid"),
        panel.grid.major.x = element_line(color = "gray", linetype = "solid"))
plot_followup_time

# Vital status
plot_vital_status <-
  ggarrange(
    ggplot(metadata, aes(vital_status)) +
      geom_bar(fill = "lightblue", col = "black") +
      labs(title = "Vital Status", x = "", y = "") +
      theme_basic,
    
    ggplot(metadata, aes(vital_status)) +
      geom_bar(fill = "lightblue", col = "black") +
      labs(subtitle = "Grouped by sex", x = "", y = "") +
      facet_wrap(~sex) +
      theme_basic,
    ncol = 1, nrow = 2)
plot_vital_status

plot_vital_status_sex <-
ggplot(metadata, aes(vital_status)) +
  geom_bar(aes(fill = sex), col = "black") +
  labs(title = "Vital Status", x = "", y = "Individuals") +
  theme_basic +
  scale_x_discrete(labels = c(
    "dead" = "Dead",
    "alive" = "Alive",
    "not reported" = "Not reported"
  )) +
  scale_y_continuous(labels = scales::comma_format()) +
  scale_fill_manual(name = "Genotypic Sex",
                    values = c(my_color_palette[4], my_color_palette[5]), # Change color (F,M)
                    labels = c("female genotypic sex" = "Female",
                               "male genotypic sex" = "Male")) +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "solid"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
plot_vital_status_sex

### Plot survival per project !!!!!!--------------------------------------------------
library(survival)
library(survminer)

metasurv <- metadata[complete.cases(metadata$days_to_death, metadata$vital_status, metadata$project), ]

# Create a survival object
survival_object <- Surv(time = metasurv$days_to_death, event = metasurv$vital_status)

# Fit a Kaplan-Meier curve
km_fit <- survfit(survival_object ~ project, data = metasurv)

# Plot the Kaplan-Meier curve
# ggsurvplot(km_fit, data = metasurv, pval = TRUE, conf.int = TRUE)










# Biosamples -------------------------------------------------------------------
## Tables ----------------------------------------------------------------------

### Table 2: Overview Samples
table2 <- metadata %>%
  select(c(substage, histological_diagnosis, icdo_morphology,
           icdo_topography, sample_origin, tumor_type, collection_moment_years))
## sample_origin, icdo_topography, icdo_morphology,
## histologicals_diagnosis too long

table2 <- metadata %>%
  select(c(substage, tumor_type, collection_moment_years, project))

CreateTableOne(data = table2)



# Substage Table
substages <- metadata %>%
  select(c(stage, substage)) %>%
  filter(stage != "Stage Unknown", stage != "Stage X",
         stage != "Stage 0", stage != "Stage I/II NOS")



### Table 3: Project
table <- metadata %>%
  select(c(sex, birth_year, age_at_diagnosis, stage,
           ethnicity, vital_status,
           days_to_death, followup_time_months, project))

table1 <- CreateTableOne(data = table)

ptable1 <- print(table1, printToggle = FALSE, noSpaces = TRUE)
kable(ptable1, format = "latex")
# table2 <- CreateTableOne(data = table, strata = c("project"))



### Histological diagnosis
##
## Maybe group into organs?
##
table3 <- metadata %>%
  select(c(histological_diagnosis, icdo_morphology, icdo_topography, sample_origin))

table3 <- CreateTableOne(data = table3)


## Visual exploration ----------------------------------------------------------

### Tumor type -----------------------------------------------------------------

levels(metadata$tumor_type) <- c("New Primary", "Metastatic", "Primary Blood",
                                 "Primary Tumor", "Recurrent Tumor")
metadata$tumor_type <- factor(
  metadata$tumor_type,
  levels = names(sort(table(metadata$tumor_type), decreasing = FALSE))
)

plot_tumor_type <-
  ggplot(metadata, aes(tumor_type)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "Tumor Type", x = "", y = "Samples (log10)") +
  scale_y_log10(labels = label_log(digits = 2)) +
  theme_basic +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "solid"))+
  coord_flip() +
  geom_text(stat = "count", aes(label = ..count..), hjust = 1.5)+
  geom_text(stat = "count", aes(label = scales::percent(round(..count../sum(..count..), 2))), hjust = -0.5, vjust = 1.5)
plot_tumor_type

### Collection moment ----------------------------------------------------------

# Make histogram with mean overlaid
plot_collection_moment <-
ggplot(metadata, aes(collection_moment_years)) +
  geom_histogram(fill = "lightblue", col = "black") +
  geom_vline(aes(xintercept = mean(collection_moment_years, na.rm = TRUE)), size = 1,
             lty = 2, col = "red", alpha = .8, show.legend = FALSE) +
  labs(
    title = "Collection Moment",
    subtitle = "Mean overlaid",
    x = "Individual Age",
    y = "Samples") +
  theme_basic +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "solid"))
plot_collection_moment

### Sample origin  -------------------------------------------------------------

# Plot the ordered bar chart
plot_frequent_origins <-
  plot_top_levels(metadata, "sample_origin", include_others = FALSE) +
  labs(title = "Most frequent sample origins", x = "", y = "Samples") +
  theme_basic +
  theme(panel.grid.major.x = element_line(color = "gray", linetype = "solid")) +
  coord_flip() +
  geom_text(stat = "count", aes(label = ..count..), hjust = 1.5)+
  geom_text(stat = "count", aes(label = scales::percent(round(..count../sum(..count..), 2))), hjust = -0.5, vjust = 1.5)
plot_frequent_origins


### Histological diagnosis -----------------------------------------------------

# Plot the top 10 histological diagnoses
plot_frequent_histologicals <-
  plot_top_levels(metadata, "histological_diagnosis",top_n = 20, include_others = F) +
  labs(title = "Most frequent Diagnoses", x = "", y = "Count")+
  theme_tight_bars

plot_frequent_histologicals

# Plot all sorted histological diagnoses
metadata$histological_diagnosis <- factor(
  metadata$histological_diagnosis,
  levels = names(sort(table(metadata$histological_diagnosis), decreasing = F))
)

ggplot(metadata, aes(histological_diagnosis)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "Histological Diagnoses", x = "", y = "Count") +
  coord_flip() +
  theme_tight_bars


### ICDO morphology ------------------------------------------------------------
# Plot the top 10 ICDO morphologies
metadata$icdo_morphology <- factor(
  metadata$icdo_morphology,
  levels = names(sort(table(metadata$icdo_morphology), decreasing = TRUE))
)

# Create a table of ICDO morphology counts and sort by count in descending order
morphology_counts <- table(metadata$icdo_morphology)
sorted_morphologies <- names(sort(morphology_counts, decreasing = TRUE))

# Identify the top 10 morphologies
top_morphologies <- sorted_morphologies[1:10]

# Create a new factor variable for ICDO morphologies, combining top morphologies and "Others"
metadata$top_morphology <- factor(
  metadata$icdo_morphology,
  levels = c(top_morphologies, "Others")
)

# Assign the "Others" label to morphologies not in the top 10
metadata$top_morphology[!(metadata$top_morphology %in% top_morphologies)] <- "Others"

# Filter the data to include only the top 10 morphologies
metadata_filtered <- metadata[metadata$icdo_morphology %in% top_morphologies, ]

# Create a new factor variable for morphologies with levels as top morphologies
metadata_filtered$icdo_morphology <- factor(
  metadata_filtered$icdo_morphology,
  levels = top_morphologies
)

# Plot the ordered bar chart
plot_frequent_morphologies <-
  ggplot(metadata_filtered, aes(icdo_morphology)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "Most frequent ICDO morphologies", x = "", y = "Count") +
  theme_basic

plot_frequent_morphologies


### ICDO topography ------------------------------------------------------------

# Plot the top 10 ICDO topographies

metadata$icdo_topography <- factor(
  metadata$icdo_topography,
  levels = names(sort(table(metadata$icdo_topography), decreasing = TRUE))
)

# Create a table of ICDO topography counts and sort by count in descending order
topography_counts <- table(metadata$icdo_topography)
sorted_topographies <- names(sort(topography_counts, decreasing = TRUE))

# Identify the top 10 topographies
top_topographies <- sorted_topographies[1:10]

# Create a new factor variable for ICDO topographies, combining top topographies and "Others"
metadata$top_topography <- factor(
  metadata$icdo_topography,
  levels = c(top_topographies, "Others")
)

# Assign the "Others" label to topographies not in the top 10
metadata$top_topography[!(metadata$top_topography %in% top_topographies)] <- "Others"

# Filter the data to include only the top 10 topographies
metadata_filtered <- metadata[metadata$icdo_topography %in% top_topographies, ]

# Create a new factor variable for topographies with levels as top topographies
metadata_filtered$icdo_topography <- factor(
  metadata_filtered$icdo_topography,
  levels = top_topographies
)

# Plot the ordered bar chart
plot_frequent_topographies <-
  ggplot(metadata_filtered, aes(icdo_topography)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "Most frequent ICDO topographies", x = "", y = "Count") +
  theme_basic

plot_frequent_topographies


### Sample origin == ICDO topography -------------------------------------------
ggarrange(plot_frequent_topographies,
          plot_frequent_origins, ncol=2)


# Results ----------------------------------------------------------------------

individuals_summary <- 
ggarrange(
  plot_sex,
  plot_stage,
  plot_ethnicities,
  plot_followup_time,
  plot_age_by_sex,
  plot_vital_status,
  ncol = 3, nrow = 2
)

sample_summary <-
ggarrange(
  plot_collection_moment,
  plot_tumor_type,
  plot_frequent_origins,
  plot_frequent_histologicals,
  plot_frequent_topographies,
  plot_frequent_morphologies,
  ncol = 3, nrow = 2
)



