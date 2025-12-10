library(tidyverse)
library(plyr)
library(dplyr)
library(readxl)
library(readr)
library(purrr)
library(ggplot2)
library(phyloseq)

#### subset data for plots ####
# load data
physeq <- readRDS("ITS.physeq_filtered.2022.RDS")

# make binary
physeq_binary <- physeq

# Apply the function to each element of the OTU table
physeq_binary <- transform_sample_counts(physeq_binary, function(x) ifelse(x>0, 1, 0))

# we want only Gala
physeq_binary <- subset_samples(physeq_binary, Cultivar == "Gala")

# make subsets for each tisue type and timepoint
physeq_pulp_fb <- subset_samples(physeq_binary, Tissue.types == "pulp" & Timepoints == "fb")
physeq_pulp_pf <- subset_samples(physeq_binary, Tissue.types == "pulp" & Timepoints == "pf")
physeq_pulp_fl1 <- subset_samples(physeq_binary, Tissue.types == "pulp" & Timepoints == "fl1")
physeq_pulp_fl2 <- subset_samples(physeq_binary, Tissue.types == "pulp" & Timepoints == "fl2")
physeq_pulp_fd2 <- subset_samples(physeq_binary, Tissue.types == "pulp" & Timepoints == "fd2")
physeq_pulp_fm1 <- subset_samples(physeq_binary, Tissue.types == "pulp" & Timepoints == "fm1")
physeq_pulp_fm2 <- subset_samples(physeq_binary, Tissue.types == "pulp" & Timepoints == "fm2")

physeq_skin_fb <- subset_samples(physeq_binary, Tissue.types == "skin" & Timepoints == "fb")
physeq_skin_pf <- subset_samples(physeq_binary, Tissue.types == "skin" & Timepoints == "pf")
physeq_skin_fl1 <- subset_samples(physeq_binary, Tissue.types == "skin" & Timepoints == "fl1")
physeq_skin_fl2 <- subset_samples(physeq_binary, Tissue.types == "skin" & Timepoints == "fl2")
physeq_skin_fd2 <- subset_samples(physeq_binary, Tissue.types == "skin" & Timepoints == "fd2")
physeq_skin_fm1 <- subset_samples(physeq_binary, Tissue.types == "skin" & Timepoints == "fm1")
physeq_skin_fm2 <- subset_samples(physeq_binary, Tissue.types == "skin" & Timepoints == "fm2")

# keep only those taxa present in 3 or more samples
physeq_pulp_fb <- prune_taxa(taxa_sums(physeq_pulp_fb) > 2, physeq_pulp_fb)
physeq_pulp_pf <- prune_taxa(taxa_sums(physeq_pulp_pf) > 2, physeq_pulp_pf)
physeq_pulp_fl1 <- prune_taxa(taxa_sums(physeq_pulp_fl1) > 2, physeq_pulp_fl1)
physeq_pulp_fl2 <- prune_taxa(taxa_sums(physeq_pulp_fl2) > 2, physeq_pulp_fl2)
physeq_pulp_fd2 <- prune_taxa(taxa_sums(physeq_pulp_fd2) > 2, physeq_pulp_fd2)
physeq_pulp_fm1 <- prune_taxa(taxa_sums(physeq_pulp_fm1) > 2, physeq_pulp_fm1)
physeq_pulp_fm2 <- prune_taxa(taxa_sums(physeq_pulp_fm2) > 2, physeq_pulp_fm2)

physeq_skin_fb <- prune_taxa(taxa_sums(physeq_skin_fb) > 2, physeq_skin_fb)
physeq_skin_pf <- prune_taxa(taxa_sums(physeq_skin_pf) > 2, physeq_skin_pf)
physeq_skin_fl1 <- prune_taxa(taxa_sums(physeq_skin_fl1) > 2, physeq_skin_fl1)
physeq_skin_fl2 <- prune_taxa(taxa_sums(physeq_skin_fl2) > 2, physeq_skin_fl2)
physeq_skin_fd2 <- prune_taxa(taxa_sums(physeq_skin_fd2) > 2, physeq_skin_fd2)
physeq_skin_fm1 <- prune_taxa(taxa_sums(physeq_skin_fm1) > 2, physeq_skin_fm1)
physeq_skin_fm2 <- prune_taxa(taxa_sums(physeq_skin_fm2) > 2, physeq_skin_fm2)

# make list of taxa present in each sample type
otu_pulp_fb <- rownames(as.data.frame(otu_table(physeq_pulp_fb)))
otu_pulp_pf <- rownames(as.data.frame(otu_table(physeq_pulp_pf)))
otu_pulp_fl1 <- rownames(as.data.frame(otu_table(physeq_pulp_fl1)))
otu_pulp_fl2 <- rownames(as.data.frame(otu_table(physeq_pulp_fl2)))
otu_pulp_fd2 <- rownames(as.data.frame(otu_table(physeq_pulp_fd2)))
otu_pulp_fm1 <- rownames(as.data.frame(otu_table(physeq_pulp_fm1)))
otu_pulp_fm2 <- rownames(as.data.frame(otu_table(physeq_pulp_fm2)))

otu_skin_fb <- rownames(as.data.frame(otu_table(physeq_skin_fb)))
otu_skin_pf <- rownames(as.data.frame(otu_table(physeq_skin_pf)))
otu_skin_fl1 <- rownames(as.data.frame(otu_table(physeq_skin_fl1)))
otu_skin_fl2 <- rownames(as.data.frame(otu_table(physeq_skin_fl2)))
otu_skin_fd2 <- rownames(as.data.frame(otu_table(physeq_skin_fd2)))
otu_skin_fm1 <- rownames(as.data.frame(otu_table(physeq_skin_fm1)))
otu_skin_fm2 <- rownames(as.data.frame(otu_table(physeq_skin_fm2)))

# make lists of overlapping taxa
pulp_fb_and_fm1 <- intersect(otu_pulp_fb, otu_pulp_fm1)
skin_fb_and_fm1 <- intersect(otu_skin_fb, otu_skin_fm1)

pulp_fb_and_fm2 <- intersect(otu_pulp_fb, otu_pulp_fm2)
skin_fb_and_fm2 <- intersect(otu_skin_fb, otu_skin_fm2)

shared_all_pulp <- Reduce(intersect, list(otu_pulp_fb,
                                          otu_pulp_pf,
                                          otu_pulp_fl1,
                                          otu_pulp_fl2,
                                          otu_pulp_fd2,
                                          otu_pulp_fm1,
                                          otu_pulp_fm2))

shared_all_skin <- Reduce(intersect, list(otu_skin_fb,
                                          otu_skin_pf,
                                          otu_skin_fl1,
                                          otu_skin_fl2,
                                          otu_skin_fd2,
                                          otu_skin_fm1,
                                          otu_skin_fm2))


#### code for function to plot ####
plot_selected_otus <- function(data, title, title2, otu_subset, physeq, label_separater, file_destination) {
  
  # calculate relative abundance
  rel_data <- data
  rel_data[2:41] <- lapply(data[2:41], function(x) (x/sum(x))*100)

  # average sample replicates so that there is one column per timepoint
  avg_data <- rel_data %>%
    dplyr::rename(otu = X) %>%  # Rename column to OTUs
    gather(key = "sample", value = "value", -1) %>%     # Reshape the data to long format
    separate(sample, into = c("sample", "replicate"), sep = label_separater) %>%  # Separate Sample and Replicate
    dplyr::group_by(sample, otu) %>%
    dplyr::summarize(
      AvgValue = mean(value, na.rm = TRUE), 
      std_dev = sd(value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop") 

  # filter based on only those OTUs of interest
  avg_data <- filter(avg_data, avg_data$otu %in% otu_subset)

  # add genus names for graph
  tax <- as.data.frame(tax_table(physeq))
  tax <- filter(tax, rownames(tax) %in% otu_subset)
  tax$otu <- rownames(tax)

  avg_data <- left_join(avg_data, tax, by = "otu")

  # visualize selected OTUs
  plot_data <- avg_data
  plot_data$sample <- factor(x = plot_data$sample, levels = c("fb", "pf", "fl1", "fl2", "fd2",  "fm1", "fm2"))
  plot_data[is.na(plot_data)] <- "Unknown"
  plot_data$otu_genus <- paste0(plot_data$Genus, " - ", plot_data$otu)

  plot_otu <- ggplot(data = plot_data, aes(x =  sample , y = AvgValue , 
                                           group = otu, col = otu_genus)) +   
    geom_point() +  
    geom_line(alpha = 1) + 
    geom_errorbar(aes(ymin = AvgValue - std_dev, ymax = AvgValue + std_dev), width = 0.2) +
    theme_bw() +  
    theme(axis.text.x = element_text(angle = 90 , vjust = 0.4)) +
    labs(x = "Timepoint", y = "Relative Abundance", title = title)
  
  plot_otu_wrap <- ggplot(data = plot_data, aes(x =  sample , y = AvgValue , 
                                           group = otu, col = otu_genus)) +   
    geom_point() +  
    geom_line(alpha = 1) + 
    geom_errorbar(aes(ymin = AvgValue - std_dev, ymax = AvgValue + std_dev), width = 0.2) +
    theme_bw() +  
    theme(axis.text.x = element_text(angle = 90 , vjust = 0.4),
          strip.text = element_text(size = 8, angle = 0)) +
    labs(x = "Timepoint", y = "Relative Abundance", title = title) + 
    facet_wrap(~otu_genus, scales = "free_y", ncol = 1)

# plot_otu_genus <- ggplot(data = plot_data, aes(x =  sample , y = AvgValue , 
#                                               group = otu, col = otu, linetype = Genus)) +   
#   geom_point() +  
#   geom_line(alpha = 1) + 
#   geom_errorbar(aes(ymin = AvgValue - std_dev, ymax = AvgValue + std_dev), width = 0.2) +
#   theme_bw() +  
#   theme(axis.text.x = element_text(angle = 90 , vjust = 0.4)) +
#   labs(x = "Timepoint", y = "Relative Abundance", title = title)

  pdf(file = file.path(file_destination, paste0(title, ".pdf")), height = 7, width = 12)
  print(plot_otu)
  dev.off()
  
  pdf(file = file.path(file_destination, paste0(title2, ".pdf")), height = 10, width = 10)
  print(plot_otu_wrap)
  dev.off()

}


#### make plots for selected OTUs ####
# ITS pulp - shared all timepoints
data <- read.csv("DESeq2_Results_ITS/pulp/normalized_data_pulp.csv")
title <- "ITS_Gala_Pulp_OTUs_Shared_all_timepoints"
title2 <- "ITS_Gala_Pulp_OTUs_Shared_all_timepoints_separated"
otu_subset <- shared_all_pulp
physeq <- readRDS("ITS.physeq_filtered.2022.rds")
label_separater <- "\\.g\\.p"
file_destination <- "SelectedOTU_Plots_ITS"

plot_selected_otus(data, title, title2, otu_subset, physeq, label_separater, file_destination)

# ITS pulp - shared fb and fm2
data <- read.csv("DESeq2_Results_ITS/pulp/normalized_data_pulp.csv")
title <- "ITS_Gala_Pulp_OTUs_Shared_fb_and_fm2"
title2 <- "ITS_Gala_Pulp_OTUs_Shared_fb_and_fm2_separated"
otu_subset <- pulp_fb_and_fm2
physeq <- readRDS("ITS.physeq_filtered.2022.rds")
label_separater <- "\\.g\\.p"
file_destination <- "SelectedOTU_Plots_ITS"

plot_selected_otus(data, title, title2, otu_subset, physeq, label_separater, file_destination)

# ITS skin shared all timepoints
data <- read.csv("DESeq2_Results_ITS/skin/normalized_data_skin.csv")
title <- "ITS_Gala_Skin_OTUs_Shared_all_timepoints"
title2 <- "ITS_Gala_Skin_OTUs_Shared_all_timepoints_separated"
otu_subset <- shared_all_skin
physeq <- readRDS("ITS.physeq_filtered.2022.rds")
label_separater <- "\\.g\\.s"
file_destination <- "SelectedOTU_Plots_ITS"

plot_selected_otus(data, title, title2, otu_subset, physeq, label_separater, file_destination)

# ITS skin shared fb and fm2
data <- read.csv("DESeq2_Results_ITS/skin/normalized_data_skin.csv")
title <- "ITS_Gala_Skin_OTUs_Shared_fb_and_fm2"
title2 <- "ITS_Gala_Skin_OTUs_Shared_fb_and_fm2_separated"
otu_subset <- skin_fb_and_fm2
physeq <- readRDS("ITS.physeq_filtered.2022.rds")
label_separater <- "\\.g\\.s"
file_destination <- "SelectedOTU_Plots_ITS"

plot_selected_otus(data, title, title2, otu_subset, physeq, label_separater, file_destination)



#### make venn diagram ####
# Install and load the venn package
install.packages("venn")
library(venn)

# Create the Venn diagram
# ITS pulp
pdf("pulp_ITS_venndiagram.pdf", height = 6, width = 8)
venn(list(fb = otu_pulp_fb,
          pf = otu_pulp_pf,
          fl1 = otu_pulp_fl1,
          fl2 = otu_pulp_fl2,
          fd2 = otu_pulp_fd2,
          fm1 = otu_pulp_fm1,
          fm2 = otu_pulp_fm2),
     counts.show = TRUE, 
     zcolor = "style",
     main = "ITS_Pulp")  
dev.off()

# 16S skin
pdf("skin_ITS_venndiagram.pdf", height = 6, width = 8)
venn(list(fb = otu_skin_fb,
          pf = otu_skin_pf,
          fl1 = otu_skin_fl1,
          fl2 = otu_skin_fl2,
          fd2 = otu_skin_fd2,
          fm1 = otu_skin_fm1,
          fm2 = otu_skin_fm2),
     counts.show = TRUE, 
     zcolor = "style",
     main = "ITS_Skin")  
dev.off()


#### pulp fb vs fm2 only - the facet wrapped plot had too many ####
#### code for function to plot ####
plot_selected_otus <- function(data, title, title2, otu_subset, physeq, label_separater, file_destination) {
  
  # calculate relative abundance
  rel_data <- data
  rel_data[2:41] <- lapply(data[2:41], function(x) (x/sum(x))*100)
  
  # average sample replicates so that there is one column per timepoint
  avg_data <- rel_data %>%
    dplyr::rename(otu = X) %>%  # Rename column to OTUs
    gather(key = "sample", value = "value", -1) %>%     # Reshape the data to long format
    separate(sample, into = c("sample", "replicate"), sep = label_separater) %>%  # Separate Sample and Replicate
    dplyr::group_by(sample, otu) %>%
    dplyr::summarize(
      AvgValue = mean(value, na.rm = TRUE), 
      std_dev = sd(value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop") 
  
  # filter based on only those OTUs of interest
  avg_data <- filter(avg_data, avg_data$otu %in% otu_subset)
  
  # add genus names for graph
  tax <- as.data.frame(tax_table(physeq))
  tax <- filter(tax, rownames(tax) %in% otu_subset)
  tax$otu <- rownames(tax)
  
  avg_data <- left_join(avg_data, tax, by = "otu")
  
  # visualize selected OTUs
  plot_data <- avg_data
  plot_data$sample <- factor(x = plot_data$sample, levels = c("fb", "pf", "fl1", "fl2", "fd2",  "fm1", "fm2"))
  plot_data[is.na(plot_data)] <- "Unknown"
  plot_data$otu_genus <- paste0(plot_data$Genus, " - ", plot_data$otu)
  
  plot_otu <- ggplot(data = plot_data, aes(x =  sample , y = AvgValue , 
                                           group = otu, col = otu_genus)) +   
    geom_point() +  
    geom_line(alpha = 1) + 
    geom_errorbar(aes(ymin = AvgValue - std_dev, ymax = AvgValue + std_dev), width = 0.2) +
    theme_bw() +  
    theme(axis.text.x = element_text(angle = 90 , vjust = 0.4)) +
    labs(x = "Timepoint", y = "Relative Abundance", title = title)
  
  plot_otu_wrap <- ggplot(data = plot_data, aes(x =  sample , y = AvgValue , 
                                                group = otu, col = otu_genus)) +   
    geom_point() +  
    geom_line(alpha = 1) + 
    geom_errorbar(aes(ymin = AvgValue - std_dev, ymax = AvgValue + std_dev), width = 0.2) +
    theme_bw() +  
    theme(axis.text.x = element_text(angle = 90 , vjust = 0.4),
          strip.text = element_text(size = 8, angle = 0)) +
    labs(x = "Timepoint", y = "Relative Abundance", title = title) + 
    facet_wrap(~otu_genus, scales = "free_y", ncol = 2)
  
  # plot_otu_genus <- ggplot(data = plot_data, aes(x =  sample , y = AvgValue , 
  #                                               group = otu, col = otu, linetype = Genus)) +   
  #   geom_point() +  
  #   geom_line(alpha = 1) + 
  #   geom_errorbar(aes(ymin = AvgValue - std_dev, ymax = AvgValue + std_dev), width = 0.2) +
  #   theme_bw() +  
  #   theme(axis.text.x = element_text(angle = 90 , vjust = 0.4)) +
  #   labs(x = "Timepoint", y = "Relative Abundance", title = title)
  
  pdf(file = file.path(file_destination, paste0(title, ".pdf")), height = 7, width = 12)
  print(plot_otu)
  dev.off()
  
  pdf(file = file.path(file_destination, paste0(title2, ".pdf")), height = 12, width = 15)
  print(plot_otu_wrap)
  dev.off()
  
}


#### function for venn diagrams ####
library(ggvenn)

# two colors
plot_venn_2 <- function(named_list, color_list, filename){
  venn_diagram <- ggvenn(
    data = named_list, 
    fill_color = color_list,
    fill_alpha = 0.5,
    set_name_size = 6,
    text_size = 6,
    auto_scale = TRUE
  )
  
  venn_diagram <- venn_diagram +
    geom_polygon(
      data = venn_diagram$data, 
      aes(x = x, y = y, color = group),
      size = 2,
      fill = NA
    ) +
    scale_color_manual(values = color_list)
  
  pdf(filename, height = 5, width = 5)
  plot(venn_diagram)
  dev.off()
}

# 4 colors
plot_venn_4 <- function(named_list, color_list, filename){
  venn_diagram <- ggvenn(
    data = named_list, 
    fill_color = color_list,
    fill_alpha = 0.5,
    set_name_size = 4,
    text_size = 3
  )
  
  venn_diagram <- venn_diagram +
    geom_polygon(
      data = venn_diagram$data, 
      aes(x = x, y = y, color = group),
      size = 1,
      fill = NA
    ) +
    scale_color_manual(values = color_list)
  
  pdf(filename, height = 5, width = 5)
  plot(venn_diagram)
  dev.off()
}

# pulp fb vs fm2
named_list <- list(Pulp_fb = otu_pulp_fb, Pulp_fm2 = otu_pulp_fm2)
color_list <- c("#0673B3", "#F1E544")
filename <- "VennDiagrams/ITS_pulp_fb_fm2.pdf"
plot_venn_2(named_list, color_list, filename)

# pulp fb, fl2, fd2, fm2
named_list <- list(Pulp_fb = otu_pulp_fb, Pulp_fl2 = otu_pulp_fl2,
                   Pulp_fd2 = otu_pulp_fd2, Pulp_fm2 = otu_pulp_fm2)
color_list <- c("#0673B3", "#F1E544", "#989898", "#D46127")
filename <- "VennDiagrams/ITS_pulp_fb__fl2_fd2_fm2.pdf"
plot_venn_4(named_list, color_list, filename)

# skin fb vs fm2
named_list <- list(Skin_fb = otu_skin_fb, Skin_fm2 = otu_skin_fm2)
color_list <- c("#0673B3", "#F1E544")
filename <- "VennDiagrams/ITS_skin_fb_fm2.pdf"
plot_venn_2(named_list, color_list, filename)

# skin fb, fl2, fd2, fm2
named_list <- list(Skin_fb = otu_skin_fb, Skin_fl2 = otu_skin_fl2,
                   Skin_fd2 = otu_skin_fd2, Skin_fm2 = otu_skin_fm2)
color_list <- c("#0673B3", "#F1E544", "#989898", "#D46127")
filename <- "VennDiagrams/ITS_skin_fb__fl2_fd2_fm2.pdf"
plot_venn_4(named_list, color_list, filename)


