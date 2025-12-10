# Load necessary libraries and functions
library(phyloseq)
library(indicspecies)


# Indicator Species Analysis function
indicator_species_analysis <- function(input_file, subset_name, output_dir) {
  # Load data
  #physeq <- readRDS(physeq)
  norm_counts <- read.csv(file = input_file, header = TRUE)
  
  # transform sample counts to relative abundance
  # calculate relative abundance
  rel_data <- norm_counts
  rel_data[2:41] <- lapply(norm_counts[2:41], function(x) (x/sum(x))*100)
  
  # transpose to put in proper format for indicspecies
  norm_counts <- rel_data
  rownames(norm_counts) <- norm_counts[,1]
  norm_counts <- norm_counts[,-1]
  norm_counts <- t(norm_counts)
  norm_counts <- as.data.frame(norm_counts)
  
  # add grouping variable for timepoint
  norm_counts <- norm_counts %>%
    mutate(Timepoint = case_when(grepl("^fb\\.", rownames(norm_counts)) ~ "fb",
                                 grepl("^pf\\.", rownames(norm_counts)) ~ "pf",
                                 grepl("^fl1\\.", rownames(norm_counts)) ~ "fl1",
                                 grepl("^fl2\\.", rownames(norm_counts)) ~ "fl2",
                                 grepl("^fd2\\.", rownames(norm_counts)) ~ "fd2",
                                 grepl("^fm1\\.", rownames(norm_counts)) ~ "fm1",
                                 grepl("^fm2\\.", rownames(norm_counts)) ~ "fm2",
                                 TRUE ~ NA_character_)  # For any other cases
    )
  
  # this code gives you indicator species for all possible groups, but not combinations of 
  # groups because there are too many and the code takes too long to run
  abund <- norm_counts[,1:(ncol(norm_counts)-1)]
  groups <- norm_counts$Timepoint
  
  inv = multipatt(abund, groups, func = "r.g", duleg = TRUE, control = how(nperm=9999))
  write.csv(capture.output(summary(inv)), file = file.path(output_dir, paste0("indsps_", subset_name, ".csv")))
}


#### Call the function for each subset ####
# pulp
input_file <- "DESeq2_Results_ITS/pulp/normalized_data_pulp.csv"
subset_name <- "gala_pulp"
output_dir <- "IndicatorSpecies_ITS"

indicator_species_analysis(input_file, subset_name, output_dir)

# skin
input_file <- "DESeq2_Results_ITS/skin/normalized_data_skin.csv"
subset_name <- "gala_skin"
output_dir <- "IndicatorSpecies_ITS"

indicator_species_analysis(input_file, subset_name, output_dir)
