library(dplyr)
library(tidyr)
library(stringr)
library(readr)

#### Format CAZy database with mapped lineages ####
# Load CAZy data
# Read the TSV file
cazy_data <- read.delim("cazy_with_lineage.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Remove JGI accessions, none of these are mapped to lineages
cazy_filtered <- subset(cazy_data, source != "jgi")

# pivot the table so that we have a count of cazy families per taxid
# First, ensure taxid is not NA (since some entries may still be missing it)
cazy_filtered_non_na <- cazy_filtered %>%
  filter(!is.na(taxid))

# Select only the columns we want to retain (taxid and taxonomy)
taxonomy_cols <- c("taxid", "organism", "domain.y", "kingdom", "phylum", "class", "order", "family", "genus", "species")

# Step 1: Get one row per taxid with taxonomy info (assuming taxonomy is consistent per taxid)
taxonomy_lookup <- cazy_filtered_non_na %>%
  select(all_of(taxonomy_cols)) %>%
  distinct()

# Group by taxid and CAZy_family and count number of annotations
cazy_summary <- cazy_filtered_non_na %>%
  group_by(taxid, CAZy_family) %>%
  summarise(count = n(), .groups = "drop")

# Pivot to wide format: one row per taxid, columns = CAZy families, values = count
cazy_wide <- cazy_summary %>%
  pivot_wider(names_from = CAZy_family, values_from = count, values_fill = 0)

# join the taxonomy back to the cazy count
cazy_wide_with_tax <- left_join(taxonomy_lookup, cazy_wide, by = "taxid") 

# replace any NAs with 0s
cazy_wide_with_tax <- cazy_wide_with_tax %>%
  mutate(across(where(is.numeric), ~replace_na(.x, 0)))

# Create a new data frame with class-level CAZy counts
# Extract column indices for each CAZy class
gh_cols <- grep("^GH", colnames(cazy_wide_with_tax), value = TRUE)
aa_cols <- grep("^AA", colnames(cazy_wide_with_tax), value = TRUE)
cbm_cols <- grep("^CBM", colnames(cazy_wide_with_tax), value = TRUE)
ce_cols <- grep("^CE", colnames(cazy_wide_with_tax), value = TRUE)
gt_cols <- grep("^GT", colnames(cazy_wide_with_tax), value = TRUE)
pl_cols <- grep("^PL", colnames(cazy_wide_with_tax), value = TRUE)

# Add summary columns to the dataframe
cazy_class_summary <- cazy_wide_with_tax %>%
  mutate(
    GH_total = rowSums(select(., all_of(gh_cols)), na.rm = TRUE),
    AA_total = rowSums(select(., all_of(aa_cols)), na.rm = TRUE),
    CBM_total = rowSums(select(., all_of(cbm_cols)), na.rm = TRUE),
    CE_total = rowSums(select(., all_of(ce_cols)), na.rm = TRUE),
    GT_total = rowSums(select(., all_of(gt_cols)), na.rm = TRUE),
    PL_total = rowSums(select(., all_of(pl_cols)), na.rm = TRUE)
  )

#### Add CAZy data to 16S and ITS data ####
# first, clean up the names in the 16S and ITS data
clean_data <- function(data) {
  # Get the data
  data <- read.csv(data)
  
  # Clean data
  data_clean <- data %>%
    mutate(
      Species = str_replace_all(Species, "_", " "),
      across(c(Kingdom, Phylum, Class, Order, Family, Genus),
             ~ str_replace(., "_1$", ""))
    )
  
  return(data_clean)
}

pulp_16S <- clean_data("CompiledData/16S_gala_pulp_compiled_data.csv")
skin_16S <- clean_data("CompiledData/16S_gala_skin_compiled_data.csv")
pulp_ITS <- clean_data("CompiledData/ITS_gala_pulp_compiled_data.csv")
skin_ITS <- clean_data("CompiledData/ITS_gala_skin_compiled_data.csv")

############ original code before writing function ###################
# add CAZy data to each condition, match by species, genus, or family
# by species: first remove NAs
# cazy_to_join_species <- filter(cazy_class_summary, !is.na(species))
# 
# pulp_16S_species <- filter(pulp_16S, !is.na(Species))
# pulp_16S_species <- filter(pulp_16S_species, Species != "uncultured bacterium")
# pulp_16S_species <- left_join(pulp_16S_species, cazy_to_join_species[,c(1:2,10,884:889)], by = c("Species" = "species"))
# 
# # sometimes an OTU matches multiple entries in the cazy table - for rows where OTU is the same, 
# # take the highest value from each CAZy summary column, concatenate taxid and organism if multiple 
# # are found for the same OTUs
# cazy_cols <- names(pulp_16S_species)[62:67]  
# 
# pulp_16S_species_simplified <- pulp_16S_species %>%
#   group_by(otu) %>%
#   summarise(
#     across(-all_of(c(cazy_cols, "taxid", "organism")), ~ first(.x)),
#     across(all_of(cazy_cols), ~ max(.x, na.rm = TRUE)),
#     taxid = paste(unique(taxid), collapse = ", "),
#     organism = paste(unique(organism), collapse = ", "),
#     .groups = "drop"
#   ) %>%
#   mutate(across(all_of(cazy_cols), ~ na_if(.x, -Inf)))
################################################################################

# function to join the cazy data to the 16S/ITS and reduce repetition
# add CAZy data to each condition, match by species, genus, or family
# by species: first remove NAs
# sometimes an OTU matches multiple entries in the cazy table - for rows where OTU is the same, 
# take the highest value from each CAZy summary column, concatenate taxid and organism if multiple 
# are found for the same OTUs
simplify_cazy_join <- function(main_df, cazy_df, main_level, cazy_level, cazy_column) {
  
  # Filter out NAs and "uncultured bacterium" from main_df
  main_df <- main_df %>%
    filter(!is.na(.data[[main_level]])) %>%
    filter(.data[[main_level]] != "uncultured bacterium")
  
  # Filter out NAs from CAZy
  cazy_df <- cazy_df %>%
    filter(!is.na(.data[[cazy_level]]))
  
  # Select only relevant columns from CAZy
  cazy_subset <- cazy_df[, c(1:2, cazy_column, 884:889)]
  
  # Join
  joined_df <- left_join(main_df, cazy_subset, by = setNames(c(cazy_level), main_level))
  
  # Identify CAZy numeric columns
  cazy_cols <- names(joined_df)[62:67]
  
  # Collapse duplicates by OTU
  simplified_df <- joined_df %>%
    group_by(otu) %>%
    summarise(
      across(-all_of(c(cazy_cols, "taxid", "organism")), ~ first(.x)),
      across(all_of(cazy_cols), ~ max(.x, na.rm = TRUE)),
      taxid = paste(unique(taxid), collapse = ", "),
      organism = paste(unique(organism), collapse = ", "),
      .groups = "drop"
    ) %>%
    mutate(across(all_of(cazy_cols), ~ na_if(.x, -Inf)))
  
  return(simplified_df)
}

# match CAZy data by species
pulp_16S_species_simplified <- simplify_cazy_join(
  main_df = pulp_16S,
  cazy_df = cazy_class_summary,
  main_level = "Species",
  cazy_level = "species",
  cazy_column = 10
  )

skin_16S_species_simplified <- simplify_cazy_join(
  main_df = skin_16S,
  cazy_df = cazy_class_summary,
  main_level = "Species",
  cazy_level = "species",
  cazy_column = 10
)

pulp_ITS <- select(pulp_ITS, -Consistent_hierarchy)
pulp_ITS_species_simplified <- simplify_cazy_join(
  main_df = pulp_ITS,
  cazy_df = cazy_class_summary,
  main_level = "Species",
  cazy_level = "species",
  cazy_column = 10
)

skin_ITS <- select(skin_ITS, -Consistent_hierarchy)
skin_ITS_species_simplified <- simplify_cazy_join(
  main_df = skin_ITS,
  cazy_df = cazy_class_summary,
  main_level = "Species",
  cazy_level = "species",
  cazy_column = 10
)

# match CAZy data by genus
pulp_16S_genus_simplified <- simplify_cazy_join(
  main_df = pulp_16S,
  cazy_df = cazy_class_summary,
  main_level = "Genus",
  cazy_level = "genus",
  cazy_column = 9
)

skin_16S_genus_simplified <- simplify_cazy_join(
  main_df = skin_16S,
  cazy_df = cazy_class_summary,
  main_level = "Genus",
  cazy_level = "genus",
  cazy_column = 9
)

# pulp_ITS <- select(pulp_ITS, -Consistent_hierarchy) # did this previously in species
pulp_ITS_genus_simplified <- simplify_cazy_join(
  main_df = pulp_ITS,
  cazy_df = cazy_class_summary,
  main_level = "Genus",
  cazy_level = "genus",
  cazy_column = 9
)

# skin_ITS <- select(skin_ITS, -Consistent_hierarchy) # did this previously in species
skin_ITS_genus_simplified <- simplify_cazy_join(
  main_df = skin_ITS,
  cazy_df = cazy_class_summary,
  main_level = "Genus",
  cazy_level = "genus",
  cazy_column = 9
)

# match CAZy data by family
pulp_16S_family_simplified <- simplify_cazy_join(
  main_df = pulp_16S,
  cazy_df = cazy_class_summary,
  main_level = "Family",
  cazy_level = "family",
  cazy_column = 8
)

skin_16S_family_simplified <- simplify_cazy_join(
  main_df = skin_16S,
  cazy_df = cazy_class_summary,
  main_level = "Family",
  cazy_level = "family",
  cazy_column = 8
)

# pulp_ITS <- select(pulp_ITS, -Consistent_hierarchy) # did this previously in species
pulp_ITS_family_simplified <- simplify_cazy_join(
  main_df = pulp_ITS,
  cazy_df = cazy_class_summary,
  main_level = "Family",
  cazy_level = "family",
  cazy_column = 8
)

# skin_ITS <- select(skin_ITS, -Consistent_hierarchy) # did this previously in species
skin_ITS_family_simplified <- simplify_cazy_join(
  main_df = skin_ITS,
  cazy_df = cazy_class_summary,
  main_level = "Family",
  cazy_level = "family",
  cazy_column = 8
)

#### write variables to files ####
# make 1 excel file for each group (i.e. skin_16S, pulp_16S, etc)
# each tab should be the annotation at the species, genus, or family level
install.packages("openxlsx")
library(openxlsx)

##### for the CAZy data with the different classes summarized #####
write.csv(cazy_class_summary, "cazy_class_summary.csv", row.names = FALSE)

##### for the simplified data frames #####
# Your simplified dataframe names
df_names <- c(
  "pulp_16S_species_simplified", 
  "pulp_16S_genus_simplified", 
  "pulp_16S_family_simplified",
  "pulp_ITS_species_simplified",
  "pulp_ITS_genus_simplified",
  "pulp_ITS_family_simplified",
  "skin_16S_species_simplified",
  "skin_16S_genus_simplified",
  "skin_16S_family_simplified",
  "skin_ITS_species_simplified",
  "skin_ITS_genus_simplified",
  "skin_ITS_family_simplified"
)

# Extract dataset prefixes, e.g. "pulp_16S"
dataset_prefixes <- unique(sub("_(species|genus|family)_simplified$", "", df_names))

# Loop through each dataset prefix
for (prefix in dataset_prefixes) {
  # Create a new workbook for this group
  wb <- createWorkbook()
  
  # Find all dataframes matching this prefix
  dfs_for_prefix <- df_names[grepl(paste0("^", prefix), df_names)]
  
  # Loop through dataframes in this group
  for (df_name in dfs_for_prefix) {
    # Extract taxonomic level for sheet name
    tax_level <- sub(paste0(prefix, "_(species|genus|family)_simplified"), "\\1", df_name)
    
    # Get dataframe object
    df <- get(df_name)
    
    # Add worksheet named by tax level
    addWorksheet(wb, sheetName = tax_level)
    
    # Write dataframe to sheet
    writeData(wb, sheet = tax_level, x = df)
  }
  
  # Save workbook with filename like "pulp_16S_simplified.xlsx"
  saveWorkbook(wb, file = paste0(prefix, "_simplified.xlsx"), overwrite = TRUE)
}


#### make plots for interpretation ####
library(dplyr)
library(ggplot2)
library(tidyr)

# # this is the original code written before writing the functions, no need to run
# # if you define the functions in the code below
# # Define the indicator groups
# indicator_groups <- list(
#   bloom = c("Indicator_fb", "Indicator_pf"),
#   fruitlet = c("Indicator_fl1", "Indicator_fl2"),
#   development = c("Indicator_fd2"),
#   maturation = c("Indicator_fm1", "Indicator_fm2")
# )
# 
# # Columns 60–65 for CAZy data
# cazy_cols <- names(pulp_16S_species_simplified)[60:65]
# 
# # Build the counts table
# counts_df <- purrr::map_dfr(names(indicator_groups), function(cat) {
#   ind_cols <- indicator_groups[[cat]]
#   
#   # Filter for "Y" in any of the indicator columns
#   has_Y <- pulp_16S_species_simplified %>%
#     filter(if_any(all_of(ind_cols), ~ .x == "Y"))
#   
#   # Count all OTUs with Y
#   total_with_Y <- nrow(has_Y)
#   
#   # Count OTUs with Y AND nonzero/non-NA in CAZy cols
#   with_cazy <- has_Y %>%
#     filter(if_any(all_of(cazy_cols), ~ !is.na(.x) & .x != 0)) %>%
#     nrow()
#   
#   tibble(
#     category = cat,
#     bar_type = c("With Y", "With Y & CAZy>0"),
#     count = c(total_with_Y, with_cazy)
#   )
# })
# 
# # Plot
# ggplot(counts_df, aes(x = category, y = count, fill = bar_type)) +
#   geom_col(position = position_dodge(width = 0.9)) +
#   geom_text(
#     aes(label = count),
#     position = position_dodge(width = 0.9),
#     vjust = -0.3,
#     size = 3
#   ) +
#   labs(title = "Number of OTUs by Fruit Development Category",
#        subtitle = "Identified to species level in 16S and CAZy database",
#        x = "Development Category", 
#        y = "Number of OTUs") +
#   scale_fill_manual(name = "Category Type",
#                     values = c("#56b4e9", "#e69f00"),
#                     labels = c("Indicator Species", "Indicator and CAZy Annotation")) +
#   theme_minimal()

#### function to plot development level by number of OTUs with indicator and CAZy annotation ####
plot_indicator_cazy <- function(df, dataset_name, taxonomic_level) {
  
  # Define the indicator groups
  indicator_groups <- list(
    bloom = c("Indicator_fb", "Indicator_pf"),
    fruitlet = c("Indicator_fl1", "Indicator_fl2"),
    development = c("Indicator_fd2"),
    maturation = c("Indicator_fm1", "Indicator_fm2")
  )
  
  # Columns 60–65 for CAZy data
  cazy_cols <- names(df)[60:65]
  
  # Build the counts table
  counts_df <- purrr::map_dfr(names(indicator_groups), function(cat) {
    ind_cols <- indicator_groups[[cat]]
    
    # Filter for "Y" in any of the indicator columns
    has_Y <- df %>%
      filter(if_any(all_of(ind_cols), ~ .x == "Y"))
    
    # Count all OTUs with Y
    total_with_Y <- nrow(has_Y)
    
    # Count OTUs with Y AND nonzero/non-NA in CAZy cols
    with_cazy <- has_Y %>%
      filter(if_any(all_of(cazy_cols), ~ !is.na(.x) & .x != 0)) %>%
      nrow()
    
    tibble(
      category = cat,
      bar_type = c("With Y", "With Y & CAZy>0"),
      count = c(total_with_Y, with_cazy)
    )
  })
  
  # Specify the order you want for the x-axis
  desired_order <- c("bloom", "fruitlet", "development", "maturation")
  
  counts_df <- counts_df %>%
    mutate(category = factor(category, levels = desired_order))
  
  # Determine ITS or 16S
  seq_type <- ifelse(grepl("ITS", dataset_name, ignore.case = TRUE), "ITS", "16S")
  
  # Plot
  plot <- ggplot(counts_df, aes(x = category, y = count, fill = bar_type)) +
    geom_col(position = position_dodge(width = 0.9)) +
    geom_text(
      aes(label = count),
      position = position_dodge(width = 0.9),
      vjust = -0.3,
      size = 3
    ) +
    labs(
      title = paste0("Number of OTUs by Fruit Development Category:", dataset_name),
      subtitle = paste("Identified to", taxonomic_level, "level in", seq_type, "and CAZy database"),
      x = "Development Category", 
      y = "Number of OTUs"
    ) +
    scale_fill_manual(
      name = "Category Type",
      values = c("#56b4e9", "#e69f00"),
      labels = c("Indicator Species", "Indicator and CAZy Annotation")
    ) +
    theme_minimal()
  
  return(plot)
}

#### function to plot k means cluster by number of OTUs with CAZy annotation ####
plot_kmeans_cazy <- function(df, dataset_name, taxonomic_level) {
  
  # Columns 60–65 for CAZy data
  cazy_cols <- names(df)[60:65]
  
  # Build counts table by KmeansCluster
  counts_df_cluster <- df %>%
    group_by(KmeansCluster) %>%
    summarise(
      total_otus = n(),
      with_cazy = sum(if_any(all_of(cazy_cols), ~ !is.na(.x) & .x != 0)),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(total_otus, with_cazy),
      names_to = "bar_type",
      values_to = "count"
    )
  
  # Determine ITS or 16S
  seq_type <- ifelse(grepl("ITS", dataset_name, ignore.case = TRUE), "ITS", "16S")
  
  # Plot
  plot <- ggplot(counts_df_cluster, aes(x = factor(KmeansCluster), y = count, fill = bar_type)) +
    geom_col(position = position_dodge(width = 0.9)) +
    geom_text(
      aes(label = count),
      position = position_dodge(width = 0.9),
      vjust = -0.3,
      size = 3
    ) +
    labs(
      title = paste0("Number of OTUs by K-means Cluster:", dataset_name),
      subtitle = paste("Identified to", taxonomic_level, "level in", seq_type, "and CAZy database"),
      x = "K-means Cluster", 
      y = "Number of OTUs"
    ) +
    scale_fill_manual(
      name = "Category Type",
      values = c("#56b4e9", "#e69f00"),
      labels = c("All OTUs", "With CAZy Annotation")
    ) +
    theme_minimal()
  
  return(plot)
}

#### function to plot indicator species and total number of cazy annotations for each developmental stage ####
plot_indicator_cazy_total <- function(df, dataset_name, taxonomic_level) {
  
  indicator_groups <- list(
    bloom = c("Indicator_fb", "Indicator_pf"),
    fruitlet = c("Indicator_fl1", "Indicator_fl2"),
    development = c("Indicator_fd2"),
    maturation = c("Indicator_fm1", "Indicator_fm2")
  )
  
  cazy_cols <- names(df)[60:65]
  
  counts_df <- purrr::map_dfr(names(indicator_groups), function(cat) {
    ind_cols <- indicator_groups[[cat]]
    
    has_Y <- df %>%
      filter(if_any(all_of(ind_cols), ~ .x == "Y"))
    
    total_cazy_annotations <- has_Y %>%
      select(all_of(cazy_cols)) %>%
      replace(is.na(.), 0) %>%
      { sum(unlist(.), na.rm = TRUE) }
    
    tibble(
      category = cat,
      total_annotations = total_cazy_annotations
    )
  })
  
  desired_order <- c("bloom", "fruitlet", "development", "maturation")
  counts_df <- counts_df %>%
    mutate(category = factor(category, levels = desired_order))
  
  seq_type <- ifelse(grepl("ITS", dataset_name, ignore.case = TRUE), "ITS", "16S")
  
  plot <- ggplot(counts_df, aes(x = category, y = total_annotations, fill = category)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = total_annotations), vjust = -0.3, size = 3) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = paste0("Total CAZy Annotations for Indicators: ", dataset_name),
      subtitle = paste("Identified to", taxonomic_level, "level in", seq_type, "and CAZy database"),
      x = "Development Category",
      y = "Total CAZy Annotations"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  return(plot)
}

#### make the plots for each condition and taxonomic level ####
# plot_indicator_cazy(df, dataset_name, taxonomic_level)
# plot_kmeans_cazy(df, dataset_name, taxonomic_level)

# List of all dataset variable names
df_names <- c(
  "pulp_16S_species_simplified", 
  "pulp_16S_genus_simplified", 
  "pulp_16S_family_simplified",
  "pulp_ITS_species_simplified",
  "pulp_ITS_genus_simplified",
  "pulp_ITS_family_simplified",
  "skin_16S_species_simplified",
  "skin_16S_genus_simplified",
  "skin_16S_family_simplified",
  "skin_ITS_species_simplified",
  "skin_ITS_genus_simplified",
  "skin_ITS_family_simplified"
)

# Map each dataframe name to its taxonomic level
tax_levels <- sub(".*_(species|genus|family)_simplified", "\\1", df_names)

# Loop through each dataframe name
for (i in seq_along(df_names)) {
  df_name <- df_names[i]
  tax_level <- tax_levels[i]
  
  # Extract dataset prefix (e.g., pulp_16S, skin_ITS) by removing taxonomic level + "_simplified"
  dataset_prefix <- sub(paste0("_", tax_level, "_simplified$"), "", df_name)
  
  # Get the actual dataframe object
  df <- get(df_name)
  
  # Create plot names
  plot_name_indicator <- paste0("plot_", sub("_simplified$", "", df_name), "_indicator")
  plot_name_kmeans <- paste0("plot_", sub("_simplified$", "", df_name), "_kmeans")
  plot_name_total <- paste0("plot_", sub("_simplified$", "", df_name), "_indicator_total_cazy")
  
  # Run plotting functions and assign plots to variables in global environment
  assign(plot_name_indicator, plot_indicator_cazy(df, dataset_prefix, tax_level), envir = .GlobalEnv)
  assign(plot_name_kmeans, plot_kmeans_cazy(df, dataset_prefix, tax_level), envir = .GlobalEnv)
  assign(plot_name_total, plot_indicator_cazy_total(df, dataset_prefix, tax_level), envir = .GlobalEnv)
}

#### save plots to pdf ####
# List all plot variable names generated in the previous step
plot_vars <- ls(pattern = "^plot_")

# Open PDF device
pdf("all_CAZy_plots.pdf", width = 10, height = 7)  # adjust size as needed

# Loop through all plot variables and print each plot to a page
for (plot_name in plot_vars) {
  print(get(plot_name))
}

# Close PDF device
dev.off()
