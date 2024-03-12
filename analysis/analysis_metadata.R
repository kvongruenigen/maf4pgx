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

# Individuals ------------------------------------------------------------------
# Filter the unique individuals
unique_individuals <- metadata[!duplicated(metadata$individual_id), ]
## Tables ----------------------------------------------------------------------
table1 <- unique_individuals %>%
  select(c(sex, birth_year, age_at_diagnosis, stage,
           ethnicity, vital_status,
           days_to_death, followup_time_months))

table1 <- CreateTableOne(data = table1)
## Why is days_to_death different from survival_time_days?
## survival_time_days is derived from days_to_death -> take days_to_death

ptable1 <- print(table1, printToggle = FALSE, noSpaces = TRUE)
kable(ptable1, format = "latex")

### Ethnicity
table(unique_individuals$ethnicity)


stat.desc(unique_individuals%>%select(c(age_at_diagnosis, followup_time_months,
                            days_to_death,collection_moment_years)))

## Visual exploration ----------------------------------------------------------

### Sex ------------------------------------------------------------------------
table(unique_individuals$sex)

plot_sex <- 
  ggplot(unique_individuals, aes(sex))+
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "Genotypic Sex", x = "", y = "Individuals") +
  scale_x_discrete(labels = c(
    "female genotypic sex" = "Female",
    "male genotypic sex" = "Male"
  )) +
  theme_basic +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "solid")) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -1.3)+
  geom_text(stat = "count", aes(label = scales::percent(round(..count../sum(..count..), 4))), vjust = 1.5)
plot_sex

# Combined with vital status
ggplot(unique_individuals, aes(sex, fill= vital_status))+
  geom_bar(aes(fill = vital_status), col = "black") +
  labs(title = "Genotypic Sex", x = "", y = "Individuals") +
  theme_basic +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "solid")) +
  scale_x_discrete(labels = c(
    "female genotypic sex" = "Female",
    "male genotypic sex" = "Male")) +
  scale_fill_manual(name = "Vital Status",
                    values = my_color_palette[2:4], # Change color (F,M)
                    labels = c("alive" = "Alive",
                               "dead" = "Dead",
                               "not reported" = "Not reported")) +  # Change legend title and colors
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5)+
  geom_text(stat = "count", aes(label = scales::percent(round(..count../sum(..count..), 2))), vjust = 1.5)

table(unique_individuals$sex)
table(unique_individuals$vital_status)
table(unique_individuals$sex, unique_individuals$vital_status)

### Stage !! -------------------------------------------------------------------

# Stage Bar Chart
plot_stage <-
  ggplot(unique_individuals, aes(stage)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "Stage", x = "", y = "Individuals") +
  theme_basic +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "solid")) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5)+
  geom_text(stat = "count", aes(label = scales::percent(round(..count../sum(..count..), 2))), vjust = 1.5)
plot_stage

# Substage Bar Chart --- Don't like the legend
ggplot(unique_individuals, aes(stage)) +
  geom_bar(aes(fill = substage), col = "black", na.rm = TRUE) +
  labs(title = "Substage", x = "", y = "Count") +
  theme_basic

# Make a table for stage per project
table_stage <- unique_individuals %>%
  select(c(stage, project))
# Create a table including NA values
table_result <- table(table_stage, useNA = "always", dnn = c("Stage", "Project"))
# Print the table
print(table_result)

# Stage 0 in 7 SKCM samples
# Stage I/II NOS in 14 SKCM samples
# Stage X occurs in 11 BRCA samples
# All CESC, DLBC, GBM, LAML, LGG, PCPG, PRAD, SARC, THYM, UCEC, UCS, OV samples are NA
# ESCA & KIRC 24-27 NA (10-13%)
# Stag IS in 49 (36%) of TGCT samples
# Stage IIC in 64 SKCM samples and 1 COAD sample
# Stage IVC in 6 THCA samples and 1 HNSC
# Stage IVB in 11 HNSC samples, 6 THCA, 3 CHOL, 2 COAD


table_stage <- unique_individuals %>%
  select(substage, project)
# Create a contingency table
table_result <- table(table_stage, useNA = "always", dnn = c("Stage", "Project"))
# Calculate the percentage of samples per project in each stage
percentage_table <- prop.table(table_result, margin = 2) * 100
rounded_percentage_table <- round(percentage_table, 0)
print(table_result) # Print the percentage table



# New plot with excluded Stage X, 0, I/II NOS and NA types
filtered_unique_individuals <- unique_individuals %>%
  filter(!stage %in% c("Stage X", "Stage 0", "Stage I/II NOS")) %>%
  group_by(project) %>%
  filter(!all(is.na(stage)))


# Create a bar chart
plot_stage <-
  ggplot(filtered_unique_individuals, aes(stage)) +
  geom_bar(fill = "lightblue", col = "black") +
  geom_bar(aes(fill = substage), col = "black", na.rm = TRUE) +
  labs(title = "Stage", x = "", y = "Samples") +
  theme_basic +
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "solid")) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5)+
  geom_text(stat = "count", aes(label = scales::percent(round(..count../sum(..count..), 2))), vjust = 1.5)
plot_stage

### Age distribution -----------------------------------------------------------

plot_age_simple <-
  ggplot(unique_individuals, aes(x = age_at_diagnosis)) +
  geom_histogram(aes(y = ..density..),binwidth = 5, fill = "lightblue", col = "black") +
  labs(title = "Age Distribution at Diagnosis", x = "Age [Years]", y = "Individuals") +
  geom_vline(aes(xintercept = mean(age_at_diagnosis, na.rm = TRUE)), size = 1,
             lty = 2, col = 'red', alpha = .8, show.legend = FALSE)+
  stat_function(fun= dnorm, #function for normal distribution
                args= list(mean(unique_individuals$age_at_diagnosis, na.rm = TRUE),# R needs mean
                           sd(unique_individuals$age_at_diagnosis, na.rm = TRUE)), # and SD
                col= "blue") +
  theme_basic+
  theme(panel.grid.major.y = element_line(color = "gray", linetype = "solid"))
plot_age_simple
# Slight left skew
skewness(unique_individuals$age_at_diagnosis, na.rm = TRUE)

# Plot age distribution per sex
plot_age_by_sex <-
  ggarrange(
    ggplot(unique_individuals, aes(age_at_diagnosis)) +
      geom_histogram(aes(y =after_stat(density)), fill = 'lightblue',
                     col = 'black',binwidth = 5)+
      geom_vline(aes(xintercept = mean(age_at_diagnosis, na.rm = TRUE)),
                 linewidth = 1, lty = 2, col = 'red', alpha = .8,
                 show.legend = FALSE) +
      labs(title = 'Age distribution',
           subtitle = 'Mean overlaid') +
      xlab('Age [years]') +
      theme_basic,
    
    ggplot(unique_individuals, aes(age_at_diagnosis))+
      geom_histogram(aes(y = after_stat(density)), fill = 'lightblue', col = 'black', binwidth = 5) +
      facet_wrap(~sex,
                 labeller = labeller(sex = c(
                   "female genotypic sex" = "Female",
                   "male genotypic sex" = "Male"))) +
      geom_vline(aes(xintercept = mean(age_at_diagnosis, na.rm = TRUE)),
                 linewidth = 1, lty = 2, col = 'red', alpha = .8,
                 show.legend = FALSE) +
      labs(title = '', subtitle = 'Grouped by sex') + 
      xlab('Age [years]') +
      theme_basic,
    ncol = 1, nrow = 2)
plot_age_by_sex

# Overall state of age at diagnosis
stats_age <- unique_individuals %>%
  summarise(mean_age = mean(age_at_diagnosis, na.rm = TRUE),
            median_age = median(age_at_diagnosis, na.rm = TRUE),
            skew = skewness(age_at_diagnosis, na.rm = TRUE),
            sd_age = sd(age_at_diagnosis, na.rm = TRUE),
            se_age = sd_age / sqrt(length(age_at_diagnosis)))


# Stats by sex
stats_age_by_sex <- unique_individuals %>%
  group_by(sex) %>%
  summarise(mean_age = mean(age_at_diagnosis, na.rm = TRUE),
            skew = skewness(age_at_diagnosis, na.rm = TRUE),
            sd_age = sd(age_at_diagnosis, na.rm = TRUE),
            se_age = sd_age / sqrt(length(age_at_diagnosis)))


# Plot age distribution per project
plot_age_per_project <-
ggplot(unique_individuals, aes(age_at_diagnosis)) +
  geom_histogram(aes(y =after_stat(density)), fill = "lightblue", col = "black",
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
  stat_function(fun= dnorm, #function for normal distribution
                args= list(mean(unique_individuals$age_at_diagnosis, na.rm = TRUE),# R needs mean
                           sd(unique_individuals$age_at_diagnosis, na.rm = TRUE)), # and SD
                n= 1e2,
                col= "blue") +
  scale_fill_manual(name = "Genotypic Sex",
                    values = c(my_color_palette[4], my_color_palette[5]), # Change color (F,M)
                    labels = c("female genotypic sex" = "Female",
                               "male genotypic sex" = "Male")) +  # Change legend title and colors
  scale_y_continuous(labels = scales::number_format(scale = 1, accuracy = 0.01))
plot_age_per_project



# Filter the unique_individuals for the interesting projects

# Identify projects with significant deviation
project_stats <- unique_individuals %>%
  group_by(project) %>%
  summarize(project_mean = mean(age_at_diagnosis, na.rm = TRUE),
            project_sd = sd(age_at_diagnosis, na.rm = TRUE),
            skew = skewness(age_at_diagnosis, na.rm = TRUE)) %>%
  arrange(desc(skew))

interesting_projects <- c("ACC", "BLCA", "CESC", "KICH", "LGG", "PCPG", "TCGT", "THCA",
                          "BLCA", "CHOL", "COAD", "ESCA", "GBM", "HNSC", "KIRP",
                          "LIHC", "LUAD", "LUSC", "MESO", "READ", "STAD", "UCEC",
                          "UCS", "PRAD")

most_skewd <- project_stats %>%
  filter(skew > 0.4 | skew < -0.4)

filtered_unique_individuals <- unique_individuals %>%
  # filter(project %in% most_skewd$project)
  filter(project %in% interesting_projects)

ggplot(filtered_unique_individuals, aes(age_at_diagnosis)) +
  geom_histogram(aes(y =after_stat(density)), fill = "lightblue", col = "black",
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
  stat_function(fun= dnorm, #function for normal distribution
                args= list(mean(unique_individuals$age_at_diagnosis, na.rm = TRUE),# R needs mean
                           sd(unique_individuals$age_at_diagnosis, na.rm = TRUE)), # and SD
                n= 1e2,
                col= "blue") +
  scale_fill_manual(name = "Genotypic Sex",
                    values = c(my_color_palette[4], my_color_palette[5]), # Change color (F,M)
                    labels = c("female genotypic sex" = "Female",
                               "male genotypic sex" = "Male")) +  # Change legend title and colors
  scale_y_continuous(labels = scales::number_format(scale = 1, accuracy = 0.01))

# Make a boxplot for age_at_diagnosis per project

ggplot(unique_individuals, aes(x = project, y = age_at_diagnosis)) +
  geom_boxplot(fill = "lightblue", col = "black") +
  labs(title = "Age at Diagnosis per Project", x = "", y = "Age [Years]") +
  geom_hline(aes(yintercept = mean(age_at_diagnosis, na.rm = TRUE)), size = 1,
             lty = 2, col = 'red', alpha = .8, show.legend = FALSE)+
  theme_basic +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::number_format(scale = 1, accuracy = 0.01))

# Assuming overall_data is your overall dataset and subpopulation_data is the subset
t_test_result <- t.test(project_stats$project_mean, mu = mean(unique_individuals$age_at_diagnosis, na.rm = TRUE), alternative = "two.sided")

# Access the p-value
p_value <- t_test_result$p.value

# Create empty vectors to store t-values and p-values
t_values <- numeric()
p_values <- numeric()
overall_mean <- mean(unique_individuals$age_at_diagnosis, na.rm = TRUE)
# Iterate through each project
for (i in seq_along(project_stats$project)) {
  # Extract project-specific data
  project_data <- unique_individuals[unique_individuals$project == project_stats$project[i], "age_at_diagnosis"]
  
  # Perform t-test
  t_test_result <- t.test(project_data, mu = overall_mean, alternative = "two.sided")
  
  # Store t-value and p-value
  t_values[i] <- t_test_result$statistic
  p_values[i] <- t_test_result$p.value
}

# Create a new data frame with project names, t-values, and p-values
result_df <- data.frame(
  project = project_stats$project,
  t_value = t_values,
  p_value = p_values
) %>%
  # round the p-values
  mutate(p_value = round(p_value, 4)) %>%
  arrange(p_value)

# Print or further analyze the result_df
print(result_df)

deviating_projects <- result_df %>%
  filter(p_value < 0.05) %>%
  arrange(project)

# Turn into simple list
deviating_projects_list <- as.character(unique(deviating_projects$project))

# Check if deviating projects are in interesting projects
deviating_projects_list %in% interesting_projects

# Show the ones that are not in interesting projects
deviating_projects_list[!deviating_projects_list %in% interesting_projects]

# Plot the deviating projects
unique_individuals_filtered <- unique_individuals %>%
  filter(project %in% deviating_projects_list)

ggplot(unique_individuals_filtered, aes(x = project, y = age_at_diagnosis)) +
  geom_boxplot(fill = "lightblue", col = "black") +
  labs(title = "Age at Diagnosis per Project", x = "", y = "Age [Years]") +
  geom_hline(aes(yintercept = mean(age_at_diagnosis, na.rm = TRUE)), size = 1,
             lty = 2, col = 'red', alpha = .8, show.legend = FALSE)+
  scale_x_discrete(limits = rev(levels(unique_individuals_filtered$project))) +
  theme_thesis +
  coord_flip() +
  scale_y_continuous(labels = scales::number_format(scale = 1, accuracy = 0.01))

# Calculate the mean age of individuals per project
mean_age_per_project <- unique_individuals %>%
  group_by(project) %>%
  summarize(mean_age = mean(age_at_diagnosis, na.rm = TRUE)) %>%
  arrange(mean_age)


### Ethnicity ------------------------------------------------------------------
# Reorder the levels of race based on frequency
unique_individuals$ethnicity_sorted <- factor(
  unique_individuals$ethnicity,
  levels = names(sort(table(unique_individuals$ethnicity),
                      decreasing = TRUE)
  )
)

levels(unique_individuals$ethnicity_sorted) <- c('white', 'not reported', 'black', 'asian',
                                       'hispanic', 'pacific islander', 'native american')

# Plot ethnicity distribution
plot_ethnicities <-
  ggplot(unique_individuals, aes(x = ethnicity_sorted)) +
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
ggplot(unique_individuals, aes(ethnicity_sorted)) +
  geom_bar(aes(fill = ethnicity_sorted)) +
  labs(title = 'Ethnicity', fill = '') +
  xlab('Ethnicity') +
  ylab('Amount') +
  theme_basic +
  scale_fill_manual(values = my_color_palette) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())



### Project --------------------------------------------------------------------
# # Remove "TCGA" and "project" from project names
# unique_individuals$project <- gsub("TCGA ", "", unique_individuals$project)
# unique_individuals$project <- gsub(" project", "", unique_individuals$project)

# Sort the projects by number of individuals
individuals_project <- factor(
  unique_individuals$project,
  levels = names(sort(table(unique_individuals$project), decreasing = FALSE))
)
table(individuals_project)

# Count individuals per project
ind_per_project <- unique_individuals %>%
  group_by(project) %>%
  summarize(individuals_project = n())%>%
  # Add percentage
  mutate(percentage = round(individuals_project / sum(individuals_project) * 100, 2))%>%
  arrange(desc(individuals_project)) %>%
  # Add cumulative percentage
  mutate(cumulative_percentage = cumsum(percentage))

ind_per_project
sum(ind_per_project$percentage[1:10])

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
unique(unique_individuals$project)

### Follow-up ------------------------------------------------------------------
# Follow-up time
plot_followup_time <-
  ggplot(unique_individuals, aes(followup_time_months))+
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
    ggplot(unique_individuals, aes(vital_status)) +
      geom_bar(fill = "lightblue", col = "black") +
      labs(title = "Vital Status", x = "", y = "") +
      theme_basic +
      theme(panel.grid.major.y = element_line(color = "gray", linetype = "solid")) +
      geom_text(stat = "count", aes(label = ..count..), vjust = -0.5)+
      geom_text(stat = "count", aes(label = scales::percent(round(..count../sum(..count..), 4))), vjust = 1.5),
    
    ggplot(unique_individuals, aes(vital_status)) +
      geom_bar(fill = "lightblue", col = "black") +
      labs(subtitle = "Grouped by sex", x = "", y = "") +
      facet_wrap(~sex) +
      theme_basic,
    ncol = 1, nrow = 2)
plot_vital_status

plot_vital_status_sex <-
ggplot(unique_individuals, aes(vital_status)) +
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

# Percentage of vital status
percentage_vital_status <- prop.table(table(unique_individuals$vital_status)) * 100
rounded_percentage_vital_status <- round(percentage_vital_status, 4)
print(rounded_percentage_vital_status)


### Plot survival per project !!!!!!------------------------------------------------
# Convert relevant columns to a survival object
surv_data <- with(unique_individuals, Surv(followup_time_months, vital_status == "dead"))

# Create a survival object with the project information
surv_project <- Surv(unique_individuals$followup_time_months, unique_individuals$vital_status == "dead")

# Create a Kaplan-Meier survival curve for each project
km_survival <- survfit(surv_project ~ project, data = unique_individuals)

# Plot the survival curves
ggsurvplot(km_survival, data = unique_individuals,
           title = "Survival Plot per Project",
           legend.labs = c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC",
                           "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML",
                           "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG",
                           "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
                           "UCEC", "UCS", "UVM"),
           legend.title = "Projects", legend = "right", palette = max_colors)
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


metadata$tumor_type_sorted <- factor(
  metadata$tumor_type,
  levels = names(sort(table(metadata$tumor_type), decreasing = FALSE))
)

plot_tumor_type <-
  ggplot(metadata, aes(tumor_type_sorted)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "Tumor Type", x = "", y = "Samples (log10)") +
  scale_y_log10(labels = label_log(digits = 2), limits= c(1,20000)) +
  scale_x_discrete(labels = c(
    "Primary Blood Derived Cancer - Peripheral Blood" = "Primary Blood",
    "Additional - New Primary" = "New Primary")) +
  theme_basic +
  theme(panel.grid.major.x = element_line(color = "gray", linetype = "solid"))+
  coord_flip() +
  geom_text(stat = "count", aes(label = after_stat(count)), hjust = 1.5)+
  geom_text(stat = "count",
            aes(label = scales::percent(round(..count../sum(..count..), 4))),
            hjust = -0.5)
plot_tumor_type

metadata %>%
  group_by(tumor_type) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

metadata %>%
  filter(tumor_type == "Metastatic") %>%
  select(project, tumor_type) %>%
  group_by(project) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

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

# Calculate mean and SD
mean(metadata$collection_moment_years, na.rm = TRUE)
sd(metadata$collection_moment_years, na.rm = TRUE)

### Sample origin  -------------------------------------------------------------
### Sample origin == ICDO topography
### ICDO topography  

# Plot the top 10 ICDO topographies

metadata$icdo_topography_specific <- factor(
  metadata$icdo_topography_specific,
  levels = names(sort(table(metadata$icdo_topography_specific), decreasing = TRUE))
)

# Create a table of ICDO topography counts and sort by count in descending order
topography_counts <- table(metadata$icdo_topography_specific)
sorted_topographies <- names(sort(topography_counts, decreasing = TRUE))

# Identify the top 10 topographies
top_topographies <- sorted_topographies[1:10]

# Create a new factor variable for ICDO topographies, combining top topographies and "Others"
metadata$top_topography <- factor(
  metadata$icdo_topography_specific,
  levels = c(top_topographies, "Others")
)

# Assign the "Others" label to topographies not in the top 10
metadata$top_topography[!(metadata$top_topography %in% top_topographies)] <- "Others"

# Filter the data to include only the top 10 topographies
metadata_filtered <- metadata[metadata$icdo_topography_specific %in% top_topographies, ]

# Create a new factor variable for topographies with levels as top topographies
metadata_filtered$icdo_topography_specific <- factor(
  metadata_filtered$icdo_topography_specific,
  levels = top_topographies
)

# Plot the ordered bar chart
plot_frequent_topographies <-
  ggplot(metadata_filtered, aes(icdo_topography_specific)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "Most frequent ICDO topographies", x = "", y = "Count") +
  theme_basic

plot_frequent_topographies

### Histological diagnosis -----------------------------------------------------

# Plot the top 10 histological diagnoses
plot_frequent_histologicals <-
  plot_top_levels(metadata, "histological_diagnosis", top_n = 20, include_others = F) +
  labs(title = "Most frequent Diagnoses", x = "", y = "Count")+
  theme_tight_bars+
  coord_flip()+
  geom_text(stat = "count", aes(label = ..count..), hjust = 1.5)

plot_frequent_histologicals

# Plot all sorted histological diagnoses
metadata$histological_diagnosis_sorted <- factor(
  metadata$histological_diagnosis,
  levels = names(sort(table(metadata$histological_diagnosis), decreasing = F))
)

ggplot(metadata, aes(histological_diagnosis)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "Histological Diagnoses", x = "", y = "Count") +
  coord_flip() +
  theme_tight_bars+
  theme(axis.text.y = element_text(size = 4))

length(unique(metadata$histological_diagnosis))

# Calculate the most frequent histological diagnoses
histological_counts <- table(metadata$histological_diagnosis)
sorted_histologicals <- sort(histological_counts, decreasing = TRUE)
# Add percentage
percentage_histologicals <- round(prop.table(histological_counts) * 100, 2)
# Combine
histologicals <- data.frame(sorted_histologicals, percentage_histologicals)
top_histologicals <- sorted_histologicals[1:10]
print(top_histologicals)

# TABLE
histologicals <- metadata %>%
  group_by(histological_diagnosis) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = (Count / sum(Count)) * 100) %>%
  arrange(desc(Count))
print(histologicals)

### ICDO morphology ------------------------------------------------------------
# Plot the top 10 ICDO morphologies
metadata$icdo_morphology_sorted <- factor(
  metadata$icdo_morphology_specific,
  levels = names(sort(table(metadata$icdo_morphology_specific), decreasing = TRUE))
)

# Create a table of ICDO morphology counts and sort by count in descending order
morphology_counts <- table(metadata$icdo_morphology_specific)
sorted_morphologies <- names(sort(morphology_counts, decreasing = TRUE))

# Identify the top 10 morphologies
top_morphologies <- sorted_morphologies[1:10]

# Create a new factor variable for ICDO morphologies, combining top morphologies and "Others"
metadata$top_morphology <- factor(
  metadata$icdo_morphology_specific,
  levels = c(top_morphologies, "Others")
)

# Assign the "Others" label to morphologies not in the top 10
metadata$top_morphology[!(metadata$top_morphology %in% top_morphologies)] <- "Others"

# Filter the data to include only the top 10 morphologies
metadata_filtered <- metadata[metadata$icdo_morphology_specific %in% top_morphologies, ]

# Create a new factor variable for morphologies with levels as top morphologies
metadata_filtered$icdo_morphology_specific <- factor(
  metadata_filtered$icdo_morphology_specific,
  levels = top_morphologies
)

# Plot the ordered bar chart
plot_frequent_morphologies <-
  ggplot(metadata_filtered, aes(icdo_morphology_specific)) +
  geom_bar(fill = "lightblue", col = "black") +
  labs(title = "Most frequent ICDO morphologies", x = "", y = "Count") +
  theme_basic

plot_frequent_morphologies


plot_top_levels(metadata, "icdo_morphology_specific", top_n = 10, include_others = F)+
  theme_tight_bars+
  coord_flip()+
  geom_text(stat = "count", aes(label = ..count..), hjust = 1.5)

# Plot all sorted ICDO morphologies
metadata$icdo_morphology_sorted <- factor(
  metadata$icdo_morphology,
  levels = names(sort(table(metadata$icdo_morphology), decreasing = F))
)
ggplot(metadata, aes(icdo_morphology_sorted, fill = project)) +
  geom_bar(aes(fill = project), col = "black", position = "stack") +
  labs(title = "ICDO Morphologies", x = "", y = "Count") +
  theme_tight_bars+
  coord_flip() +
  scale_fill_manual(name = "Projects",
                    values = colors)

length(unique(metadata$icdo_morphology))

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
