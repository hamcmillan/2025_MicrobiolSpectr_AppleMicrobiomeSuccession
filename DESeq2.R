library(tidyverse)
library(plyr)
library(dplyr)
library(readxl)
library(readr)
library(purrr)
library(ggplot2)


BiocManager::install("phyloseq")
library(phyloseq)

#### load data and get rid of samples with no reads ####
# load data
physeq <- readRDS("16S.physeq_filtered.RDS")

# we want only Gala
physeq <- subset_samples(physeq, Cultivar == "Gala")

# make pulp and skin datasets
physeq_pulp <- subset_samples(physeq, Tissue.types == "pulp")
physeq_skin <- subset_samples(physeq, Tissue.types == "skin")

# remove OTUs with no reads
# prune OTUs that are not present in at least one sample
physeq_pulp <- prune_taxa(taxa_sums(physeq_pulp) > 0, physeq_pulp)
physeq_skin <- prune_taxa(taxa_sums(physeq_skin) > 0, physeq_skin)

# do any samples now have no reads/ASVs?
sum(sample_sums(physeq_pulp) == 0)
sum(sample_sums(physeq_skin) == 0)


#### PULP ####
physeq <- physeq_pulp

#### convert to DESeq2 object and run ####
library(DESeq2)

# convert any variables to factors
physeq@sam_data$Replicates_fact <- factor(physeq@sam_data$Replicates)
physeq@sam_data$Timepoints_fact <- factor(physeq@sam_data$Timepoints, ordered = FALSE, 
                                          levels = c("fb", "pf", "fl1", "fl2", "fd2",  "fm1", "fm2"))

# set reference levels
# https://www.biostars.org/p/357464/
physeq@sam_data$Replicates_fact <- relevel(physeq@sam_data$Replicates_fact, ref = 1)
#physeq@sam_data$Timepoints_fact <- relevel(physeq@sam_data$Timepoints_fact, ref = "fb") # this should be taken care of by specifying levels above

# add pseudocount to avoid zeros
physeq@otu_table <- physeq@otu_table + 1

# make DESq2 matrix and run model
# put any possible batch effects first in the model so that they are subtracted from downstream comparisons 
ddsMat <- phyloseq_to_deseq2(physeq, ~ Replicates_fact + Timepoints_fact)
dds <- DESeq(ddsMat, test="Wald", fitType="mean")

# view model information/parameters
matrix(resultsNames(dds)) #view the comparison outputs from the model

# [1,] "Intercept"                
# [2,] "Replicates_fact_2_vs_1"   
# [3,] "Replicates_fact_3_vs_1"   
# [4,] "Replicates_fact_4_vs_1"   
# [5,] "Replicates_fact_5_vs_1"   
# [6,] "Replicates_fact_6_vs_1"   
# [7,] "Timepoints_fact_pf_vs_fb" 
# [8,] "Timepoints_fact_fl1_vs_fb"
# [9,] "Timepoints_fact_fl2_vs_fb"
# [10,] "Timepoints_fact_fd2_vs_fb"
# [11,] "Timepoints_fact_fm1_vs_fb"
# [12,] "Timepoints_fact_fm2_vs_fb"

attr(dds, "modelMatrixType") #view the model type used by DESeq2

# "standard"

###### batch effects ######
# first check - is there an effect of experimental replicate that we need to block for
# there are experimental effects for some of the genes - for some downstream analysis we may 
# need to use variance stabilize the counts, apply limma's removeBatchEffect to assay(vsd), 
# then use plotPCA to plot the residuals
# https://support.bioconductor.org/p/76099/#93176
res <- results(dds, contrast = c("Replicates_fact", "2", "1")) # only 1 OTU significant
summary(res)
# out of 590 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 1, 0.17%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

res <- results(dds, contrast = c("Replicates_fact", "3", "1"))
summary(res)
# out of 590 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 19, 3.2%
# LFC < 0 (down)     : 13, 2.2%
# outliers [1]       : 0, 0%
# low counts [2]     : 389, 66%
# (mean count < 2)

res <- results(dds, contrast = c("Replicates_fact", "4", "1"))
summary(res)
# out of 590 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 6, 1%
# LFC < 0 (down)     : 7, 1.2%
# outliers [1]       : 0, 0%
# low counts [2]     : 446, 76%
# (mean count < 3)

res <- results(dds, contrast = c("Replicates_fact", "5", "1"))
summary(res)
# out of 590 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 24, 4.1%
# LFC < 0 (down)     : 18, 3.1%
# outliers [1]       : 0, 0%
# low counts [2]     : 423, 72%
# (mean count < 3)

res <- results(dds, contrast = c("Replicates_fact", "6", "1"))
summary(res)
# out of 590 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 12, 2%
# LFC < 0 (down)     : 12, 2%
# outliers [1]       : 0, 0%
# low counts [2]     : 446, 76%
# (mean count < 3)

###### main effects ######
# next, follow example 3 in ?results to get results/comparisons of interest
# timepoints effect pf vs fb
res <- results(dds, contrast=c("Timepoints_fact", "pf", "fb"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_pf_v_fb.csv")

# make df to collect results
all_Wald_results <- as.data.frame(res[,c(2,6)])
colnames(all_Wald_results) <- c("log2FoldChange_time_pf_v_fb", "padj_time_pf_v_fb")

# timepoints effect fl1 vs fb
res <- results(dds, contrast=c("Timepoints_fact", "fl1", "fb"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fl1_v_fb.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fl1_v_fb", "padj__time_fl1_v_fb")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fl2 vs fb
res <- results(dds, contrast=c("Timepoints_fact", "fl2", "fb"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fl2_v_fb.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fl2_v_fb", "padj__time_fl2_v_fb")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fd2 vs fb
res <- results(dds, contrast=c("Timepoints_fact", "fd2", "fb"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fd2_v_fb.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fd2_v_fb", "padj__time_fd2_v_fb")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm1 vs fb
res <- results(dds, contrast=c("Timepoints_fact", "fm1", "fb"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm1_v_fb.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm1_v_fb", "padj__time_fm1_v_fb")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm2 vs fb
res <- results(dds, contrast=c("Timepoints_fact", "fm2", "fb"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm2_v_fb.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm2_v_fb", "padj__time_fm2_v_fb")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fl1 vs pf
res <- results(dds, contrast=c("Timepoints_fact", "fl1", "pf"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fl1_v_pf.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fl1_v_pf", "padj__time_fl1_v_pf")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fl2 vs pf
res <- results(dds, contrast=c("Timepoints_fact", "fl2", "pf"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fl2_v_pf.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fl2_v_pf", "padj__time_fl2_v_pf")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fd2 vs pf
res <- results(dds, contrast=c("Timepoints_fact", "fd2", "pf"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fd2_v_pf.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fd2_v_pf", "padj__time_fd2_v_pf")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm1 vs pf
res <- results(dds, contrast=c("Timepoints_fact", "fm1", "pf"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm1_v_pf.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm1_v_pf", "padj__time_fm1_v_pf")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm2 vs pf
res <- results(dds, contrast=c("Timepoints_fact", "fm2", "pf"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm2_v_pf.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm2_v_pf", "padj__time_fm2_v_pf")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fl2 vs fl1
res <- results(dds, contrast=c("Timepoints_fact", "fl2", "fl1"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fl2_v_fl1.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fl2_v_fl1", "padj__time_fl2_v_fl1")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fd2 vs fl1
res <- results(dds, contrast=c("Timepoints_fact", "fd2", "fl1"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fd2_v_fl1.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fd2_v_fl1", "padj__time_fd2_v_fl1")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm1 vs fl1
res <- results(dds, contrast=c("Timepoints_fact", "fm1", "fl1"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm1_v_fl1.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm1_v_fl1", "padj__time_fm1_v_fl1")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm2 vs fl1
res <- results(dds, contrast=c("Timepoints_fact", "fm2", "fl1"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm2_v_fl1.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm2_v_fl1", "padj__time_fm2_v_fl1")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fd2 vs fl2
res <- results(dds, contrast=c("Timepoints_fact", "fd2", "fl2"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fd2_v_fl2.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fd2_v_fl2", "padj__time_fd2_v_fl2")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm1 vs fl2
res <- results(dds, contrast=c("Timepoints_fact", "fm1", "fl2"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm1_v_fl2.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm1_v_fl2", "padj__time_fm1_v_fl2")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm2 vs fl2
res <- results(dds, contrast=c("Timepoints_fact", "fm2", "fl2"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm2_v_fl2.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm2_v_fl2", "padj__time_fm2_v_fl2")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm1 vs fd2
res <- results(dds, contrast=c("Timepoints_fact", "fm1", "fd2"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm1_v_fd2.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm1_v_fd2", "padj__time_fm1_v_fd2")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm2 vs fd2
res <- results(dds, contrast=c("Timepoints_fact", "fm2", "fd2"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm2_v_fd2.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm2_v_fd2", "padj__time_fm2_v_fd2")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm2 vs fm1
res <- results(dds, contrast=c("Timepoints_fact", "fm2", "fm1"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm2_v_fm1.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm2_v_fm1", "padj__time_fm2_v_fm1")
all_Wald_results <- cbind(all_Wald_results, res)

##### Save Results #####
# we want to make sure we have a few things from this analysis:
# 1. a table with the transformed counts for each asv
# 2. a table with the asv and the padj and lfc value for each effect test
# 3. a table with the asv and associated taxonomy
# 4. a table with the normalized counts used in DESeq2 analysis

# 1
# several ways to extract transformed data, use this: https://support.bioconductor.org/p/123969/
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vsd_data <- as.data.frame(assay(vsd))
write.csv(vsd_data, file="DESeq2_Results/variance_stabilized_transform_pulp.csv")

# 2
write.csv(all_Wald_results, file="DESeq2_Results/all_Wald_results_pulp.csv")

# 3
tax <- as.data.frame(physeq@tax_table)
write.csv(tax, file="DESeq2_Results/tax_table_pulp.csv")

# 4
normalized_counts <- as.data.frame(counts(dds, normalized = T))
write.csv(normalized_counts, file="DESeq2_Results/normalized_data_pulp.csv")

# save the dds object so that it doesn't get overwritten when analyzing total
dds_pulp <- dds


#### k means clustering and plots ####
# https://www.biostars.org/p/343055/#343128
# need to start with df that has first column is geneName and 
# all other columns are relative expression in given sample
for_clust <- normalized_counts

# average sample replicates so that there is one column per timepoint
averaged_df <- for_clust %>%
  rownames_to_column(var = "otu") %>%  # Create a new column for row names (OTUs)
  gather(key = "sample", value = "value", -1) %>%     # Reshape the data to long format
  separate(sample, into = c("sample", "replicate"), sep = ".g.p") %>%  # Separate Sample and Replicate
  dplyr::group_by(sample, otu) %>%
  dplyr::summarize(AvgValue = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = AvgValue)

# calculate relative abundance
averaged_df[2:8] <- lapply(averaged_df[2:8], function(x) (x/sum(x))*100)


# kmeans
max_itr <-  50
n_clust  <-  6  ## number of cluster 
set.seed(123) ## reproduce the cluster 
kmeans_out  <- kmeans(averaged_df[2:8],n_clust,iter.max = max_itr)

# add cluster info to orig matrix 
data_with_clust_info <- averaged_df %>% 
  mutate(clust = paste("clust_", kmeans_out$cluster,sep = ""))

#data_with_clust_info <- rownames_to_column(data_with_clust_info, "otu")

# visualise  each cluster 
cluster_plot <- data_with_clust_info %>% 
  gather(key = "variable" , value = "value", -c(1,9)) %>%  ### 1 is the index of column 'otu' and 43 is the index of column 'clust'
  group_by(variable) %>%  
  dplyr::mutate(row_num =  1:n()) 
cluster_plot$variable <- factor(x = cluster_plot$variable, levels = c("fb", "pf", "fl1", "fl2", "fd2",  "fm1", "fm2"))

cluster_plot <- ggplot(data = cluster_plot, aes(x =  variable , y = value , group = row_num)) +   
  geom_point() +  
  geom_line(alpha = 1 , aes(col = as.character(clust))) + 
  theme_bw() +  
  theme(legend.position = "none" , axis.text.x = element_text(angle = 90 , vjust = 0.4)) +
  facet_wrap(~clust)

# save plot and clusters to file
write.csv(data_with_clust_info, "Figures/clusters_pulp_relabund.csv")

pdf("Figures/clusters_pulp_relabund.pdf")
cluster_plot
dev.off()


#### SKIN ####
physeq <- physeq_skin

#### convert to DESeq2 object and run ####
library(DESeq2)

# convert any variables to factors
physeq@sam_data$Replicates_fact <- factor(physeq@sam_data$Replicates)
physeq@sam_data$Timepoints_fact <- factor(physeq@sam_data$Timepoints, ordered = FALSE, 
                                          levels = c("fb", "pf", "fl1", "fl2", "fd2",  "fm1", "fm2"))

# set reference levels
# https://www.biostars.org/p/357464/
physeq@sam_data$Replicates_fact <- relevel(physeq@sam_data$Replicates_fact, ref = 1)
#physeq@sam_data$Timepoints_fact <- relevel(physeq@sam_data$Timepoints_fact, ref = "fb") # this should be taken care of by specifying levels above

# add pseudocount to avoid zeros
physeq@otu_table <- physeq@otu_table + 1

# make DESq2 matrix and run model
# put any possible batch effects first in the model so that they are subtracted from downstream comparisons 
ddsMat <- phyloseq_to_deseq2(physeq, ~ Replicates_fact + Timepoints_fact)
dds <- DESeq(ddsMat, test="Wald", fitType="mean")

# view model information/parameters
matrix(resultsNames(dds)) #view the comparison outputs from the model

# [1,] "Intercept"                
# [2,] "Replicates_fact_2_vs_1"   
# [3,] "Replicates_fact_3_vs_1"   
# [4,] "Replicates_fact_4_vs_1"   
# [5,] "Replicates_fact_5_vs_1"   
# [6,] "Replicates_fact_6_vs_1"   
# [7,] "Timepoints_fact_pf_vs_fb" 
# [8,] "Timepoints_fact_fl1_vs_fb"
# [9,] "Timepoints_fact_fl2_vs_fb"
# [10,] "Timepoints_fact_fd2_vs_fb"
# [11,] "Timepoints_fact_fm1_vs_fb"
# [12,] "Timepoints_fact_fm2_vs_fb"

attr(dds, "modelMatrixType") #view the model type used by DESeq2

# "standard"

###### batch effects ######
# first check - is there an effect of experimental replicate that we need to block for
# there are experimental effects for some of the genes - for some downstream analysis we may 
# need to use variance stabilize the counts, apply limma's removeBatchEffect to assay(vsd), 
# then use plotPCA to plot the residuals
# https://support.bioconductor.org/p/76099/#93176
res <- results(dds, contrast = c("Replicates_fact", "2", "1")) # only 1 OTU significant
summary(res)
# out of 637 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1, 0.16%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

res <- results(dds, contrast = c("Replicates_fact", "3", "1"))
summary(res)
# out of 637 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)

res <- results(dds, contrast = c("Replicates_fact", "4", "1"))
summary(res)
# out of 637 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 1, 0.16%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)

res <- results(dds, contrast = c("Replicates_fact", "5", "1"))
summary(res)
# out of 637 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 13, 2%
# LFC < 0 (down)     : 3, 0.47%
# outliers [1]       : 0, 0%
# low counts [2]     : 494, 78%
# (mean count < 2)

res <- results(dds, contrast = c("Replicates_fact", "6", "1"))
summary(res)
# out of 637 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 6, 0.94%
# LFC < 0 (down)     : 3, 0.47%
# outliers [1]       : 0, 0%
# low counts [2]     : 531, 83%
# (mean count < 3)

###### main effects ######
# next, follow example 3 in ?results to get results/comparisons of interest
# timepoints effect pf vs fb
res <- results(dds, contrast=c("Timepoints_fact", "pf", "fb"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_pf_v_fb.csv")

# make df to collect results
all_Wald_results <- as.data.frame(res[,c(2,6)])
colnames(all_Wald_results) <- c("log2FoldChange_time_pf_v_fb", "padj_time_pf_v_fb")

# timepoints effect fl1 vs fb
res <- results(dds, contrast=c("Timepoints_fact", "fl1", "fb"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fl1_v_fb.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fl1_v_fb", "padj__time_fl1_v_fb")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fl2 vs fb
res <- results(dds, contrast=c("Timepoints_fact", "fl2", "fb"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fl2_v_fb.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fl2_v_fb", "padj__time_fl2_v_fb")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fd2 vs fb
res <- results(dds, contrast=c("Timepoints_fact", "fd2", "fb"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fd2_v_fb.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fd2_v_fb", "padj__time_fd2_v_fb")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm1 vs fb
res <- results(dds, contrast=c("Timepoints_fact", "fm1", "fb"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm1_v_fb.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm1_v_fb", "padj__time_fm1_v_fb")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm2 vs fb
res <- results(dds, contrast=c("Timepoints_fact", "fm2", "fb"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm2_v_fb.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm2_v_fb", "padj__time_fm2_v_fb")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fl1 vs pf
res <- results(dds, contrast=c("Timepoints_fact", "fl1", "pf"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fl1_v_pf.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fl1_v_pf", "padj__time_fl1_v_pf")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fl2 vs pf
res <- results(dds, contrast=c("Timepoints_fact", "fl2", "pf"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fl2_v_pf.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fl2_v_pf", "padj__time_fl2_v_pf")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fd2 vs pf
res <- results(dds, contrast=c("Timepoints_fact", "fd2", "pf"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fd2_v_pf.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fd2_v_pf", "padj__time_fd2_v_pf")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm1 vs pf
res <- results(dds, contrast=c("Timepoints_fact", "fm1", "pf"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm1_v_pf.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm1_v_pf", "padj__time_fm1_v_pf")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm2 vs pf
res <- results(dds, contrast=c("Timepoints_fact", "fm2", "pf"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm2_v_pf.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm2_v_pf", "padj__time_fm2_v_pf")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fl2 vs fl1
res <- results(dds, contrast=c("Timepoints_fact", "fl2", "fl1"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fl2_v_fl1.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fl2_v_fl1", "padj__time_fl2_v_fl1")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fd2 vs fl1
res <- results(dds, contrast=c("Timepoints_fact", "fd2", "fl1"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fd2_v_fl1.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fd2_v_fl1", "padj__time_fd2_v_fl1")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm1 vs fl1
res <- results(dds, contrast=c("Timepoints_fact", "fm1", "fl1"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm1_v_fl1.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm1_v_fl1", "padj__time_fm1_v_fl1")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm2 vs fl1
res <- results(dds, contrast=c("Timepoints_fact", "fm2", "fl1"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm2_v_fl1.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm2_v_fl1", "padj__time_fm2_v_fl1")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fd2 vs fl2
res <- results(dds, contrast=c("Timepoints_fact", "fd2", "fl2"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fd2_v_fl2.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fd2_v_fl2", "padj__time_fd2_v_fl2")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm1 vs fl2
res <- results(dds, contrast=c("Timepoints_fact", "fm1", "fl2"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm1_v_fl2.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm1_v_fl2", "padj__time_fm1_v_fl2")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm2 vs fl2
res <- results(dds, contrast=c("Timepoints_fact", "fm2", "fl2"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm2_v_fl2.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm2_v_fl2", "padj__time_fm2_v_fl2")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm1 vs fd2
res <- results(dds, contrast=c("Timepoints_fact", "fm1", "fd2"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm1_v_fd2.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm1_v_fd2", "padj__time_fm1_v_fd2")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm2 vs fd2
res <- results(dds, contrast=c("Timepoints_fact", "fm2", "fd2"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm2_v_fd2.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm2_v_fd2", "padj__time_fm2_v_fd2")
all_Wald_results <- cbind(all_Wald_results, res)

# timepoints effect fm2 vs fm1
res <- results(dds, contrast=c("Timepoints_fact", "fm2", "fm1"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.csv(res, file="DESeq2_Results/time_effect_fm2_v_fm1.csv")
res <- res[,c(2,6)][match(rownames(all_Wald_results), rownames(res)),]
colnames(res) <- c("log2FoldChange_time_fm2_v_fm1", "padj__time_fm2_v_fm1")
all_Wald_results <- cbind(all_Wald_results, res)

##### Save Results #####
# we want to make sure we have a few things from this analysis:
# 1. a table with the transformed counts for each asv
# 2. a table with the asv and the padj and lfc value for each effect test
# 3. a table with the asv and associated taxonomy
# 4. a table with the normalized counts used in DESeq2 analysis

# 1
# several ways to extract transformed data, use this: https://support.bioconductor.org/p/123969/
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vsd_data <- as.data.frame(assay(vsd))
write.csv(vsd_data, file="DESeq2_Results/variance_stabilized_transform_skin.csv")

# 2
write.csv(all_Wald_results, file="DESeq2_Results/all_Wald_results_skin.csv")

# 3
tax <- as.data.frame(physeq@tax_table)
write.csv(tax, file="DESeq2_Results/tax_table_skin.csv")

# 4
normalized_counts <- as.data.frame(counts(dds, normalized = T))
write.csv(normalized_counts, file="DESeq2_Results/normalized_data_skin.csv")

# save the dds object so that it doesn't get overwritten when analyzing total
dds_skin <- dds


#### k means clustering and plots ####
# https://www.biostars.org/p/343055/#343128
# need to start with df that has first column is geneName and 
# all other columns are relative expression in given sample
for_clust <- normalized_counts

# average sample replicates so that there is one column per timepoint
averaged_df <- for_clust %>%
  rownames_to_column(var = "otu") %>%  # Create a new column for row names (OTUs)
  gather(key = "sample", value = "value", -1) %>%     # Reshape the data to long format
  separate(sample, into = c("sample", "replicate"), sep = ".g.s") %>%  # Separate Sample and Replicate
  dplyr::group_by(sample, otu) %>%
  dplyr::summarize(AvgValue = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = AvgValue)

# calculate relative abundance
averaged_df[2:8] <- lapply(averaged_df[2:8], function(x) (x/sum(x))*100)


# kmeans
max_itr <-  50
n_clust  <-  6  ## number of cluster 
set.seed(123) ## reproduce the cluster 
kmeans_out  <- kmeans(averaged_df[2:8],n_clust,iter.max = max_itr)

# add cluster info to orig matrix 
data_with_clust_info <- averaged_df %>% 
  mutate(clust = paste("clust_", kmeans_out$cluster,sep = ""))

#data_with_clust_info <- rownames_to_column(data_with_clust_info, "otu")

# visualise  each cluster 
cluster_plot <- data_with_clust_info %>% 
  gather(key = "variable" , value = "value", -c(1,9)) %>%  ### 1 is the index of column 'otu' and 43 is the index of column 'clust'
  group_by(variable) %>%  
  dplyr::mutate(row_num =  1:n()) 
cluster_plot$variable <- factor(x = cluster_plot$variable, levels = c("fb", "pf", "fl1", "fl2", "fd2",  "fm1", "fm2"))

cluster_plot <- ggplot(data = cluster_plot, aes(x =  variable , y = value , group = row_num)) +   
  geom_point() +  
  geom_line(alpha = 1 , aes(col = as.character(clust))) + 
  theme_bw() +  
  theme(legend.position = "none" , axis.text.x = element_text(angle = 90 , vjust = 0.4)) +
  facet_wrap(~clust)

# save plot and clusters to file
write.csv(data_with_clust_info, "Figures/clusters_skin_relabund.csv")

pdf("Figures/clusters_skin_relabund.pdf")
cluster_plot
dev.off()














