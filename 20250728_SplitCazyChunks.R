library(dplyr)

# Load CAZy data
cazy <- read.table("cazy_data_20250707.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                   quote = "", fill = TRUE, comment.char = "",
                   col.names = c("Family", "Domain", "Organism", "Accession", "Source"))

# Filter only NCBI entries
cazy_ncbi <- cazy %>% filter(Source == "ncbi")

# Get unique accessions
unique_ncbi_accessions <- unique(cazy_ncbi$Accession)

# Split into chunks of 10,000
chunk_size <- 10000
chunks <- split(unique_ncbi_accessions, ceiling(seq_along(unique_ncbi_accessions) / chunk_size))

dir.create("cazy_taxonomy_chunks", showWarnings = FALSE)

# Save each chunk as a plain text file
for (i in seq_along(chunks)) {
  writeLines(chunks[[i]], sprintf("cazy_taxonomy_chunks/accessions_%04d.txt", i))
}
