install.packages("rentrez")
install.packages("xml2")

library(rentrez)
library(xml2)
library(dplyr)

#### Function to check accession formatting ahead of the main function ####
# this excludes a lot of accessions that actually do have entries, so I wouldn't recommend
check_accession_format <- function(accessions) {
  # Match: 1+ uppercase letters or underscore, digits, then .version (e.g. WP_123456789.1, AAB12345.1)
  valid_pattern <- "^[A-Z_]+[0-9]+\\.[0-9]+$"
  invalid <- accessions[!grepl(valid_pattern, accessions)]
  
  if (length(invalid) > 0) {
    warning("The following accessions may be malformed or missing version numbers:\n",
            paste(invalid, collapse = ", "))
  } else {
    message("All accessions appear properly formatted.")
  }
  return(invalid)
}

#### Function to get taxonomy from NCBI ####
get_taxonomy_from_accession_full <- function(accession) {
  tryCatch({
    search_res <- entrez_search(db = "protein", term = paste0(accession, "[ACCN]"))
    if (length(search_res$ids) == 0) {
      warning(paste("No protein ID found for", accession))
      stop("Search failed")
    }
    protein_id <- search_res$ids[1]
    
    sum_res <- entrez_summary(db = "protein", id = protein_id)
    taxid <- sum_res$taxid
    organism <- sum_res$organism
    
    if (is.null(taxid)) {
      message(paste("No taxid found for", accession, "- attempting elink to nucleotide"))
      links <- entrez_link(dbfrom = "protein", id = protein_id, db = "nucleotide")
      if (!is.null(links$links$protein_nucleotide)) {
        nucl_id <- links$links$protein_nucleotide[1]
        nucl_summary <- entrez_summary(db = "nucleotide", id = nucl_id)
        taxid <- nucl_summary$taxid
        organism <- nucl_summary$organism
      } else {
        warning(paste("No linked nucleotide record found for", accession))
        stop("TaxID missing")
      }
    }
    
    tax_raw <- entrez_fetch(db = "taxonomy", id = taxid, rettype = "xml", parsed = FALSE)
    tax_xml <- read_xml(tax_raw)
    
    lineage_nodes <- xml_find_all(tax_xml, ".//LineageEx/Taxon")
    ranks <- xml_text(xml_find_all(lineage_nodes, "./Rank"))
    names <- xml_text(xml_find_all(lineage_nodes, "./ScientificName"))
    lineage <- setNames(names, ranks)
    
    current_name <- xml_text(xml_find_first(tax_xml, ".//Taxon/ScientificName"))
    current_rank <- xml_text(xml_find_first(tax_xml, ".//Taxon/Rank"))
    lineage[current_rank] <- current_name
    
    major_ranks <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
    major_lineage <- setNames(rep(NA_character_, length(major_ranks)), major_ranks)
    for (rank in major_ranks) {
      if (!is.na(lineage[rank])) {
        major_lineage[rank] <- lineage[rank]
      }
    }
    
    return(data.frame(
      accession = accession,
      taxid = taxid,
      scientific_name = organism,
      t(major_lineage),
      stringsAsFactors = FALSE,
      row.names = NULL
    ))
    
  }, error = function(e) {
    warning(paste("Failed to retrieve taxonomy for", accession, "-", e$message))
    return(data.frame(
      accession = accession,
      taxid = NA,
      scientific_name = NA,
      superkingdom = NA,
      kingdom = NA,
      phylum = NA,
      class = NA,
      order = NA,
      family = NA,
      genus = NA,
      species = NA,
      stringsAsFactors = FALSE,
      row.names = NULL
    ))
  })
}

# # Test with example accession:
# res <- get_taxonomy_from_accession_full("UBD70155.1")
# res <- get_taxonomy_from_accession_full("A43802")
# print(res)

#### Use function to add full NCBI taxonomy to the CAZy database ####
# Load CAZy data
cazy <- read.table("cazy_data_20250707.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                   quote = "", fill = TRUE, comment.char = "",
                   col.names = c("Family", "Domain", "Organism", "Accession", "Source"))

# Filter only NCBI entries
cazy_ncbi <- cazy %>% filter(Source == "ncbi")

# filter JGI entries
cazy_jgi <- cazy %>% filter(Source == "jgi")

# Get unique accessions
unique_ncbi_accessions <- unique(cazy_ncbi$Accession)
unique_jgi_accessions <- unique(cazy_jgi$Accession)

# Check accession formatting
invalid_accessions <- check_accession_format(unique_ncbi_accessions)

# Lookup taxonomy for all accessions (takes time!)
# Store failures
failed_accessions <- c()

# Loop with progress logging
taxonomy_results <- lapply(unique_ncbi_accessions, function(acc) {
  res <- get_taxonomy_from_accession_full(acc)
  if (is.na(res$taxid)) {
    failed_accessions <<- c(failed_accessions, acc)
  }
  return(res)
})

# Combine into one dataframe
taxonomy_df <- do.call(rbind, taxonomy_results)

# taxonomy_results <- lapply(unique_ncbi_accessions, get_taxonomy_from_accession_full)
# 
# # Combine into a dataframe
# taxonomy_df <- bind_rows(taxonomy_results)

# Join taxonomy back into CAZy table
cazy_ncbi_annotated <- left_join(cazy_ncbi, taxonomy_df, by = "Accession")


#### test an individual accession ####
rentrez::set_entrez_key("3efcf5827b1ed46c10b3a518d11bd70bcc09")

# Try searching for the accession in protein DB
search_res <- entrez_search(db = "protein", term = "ABM99806.1[ACCN]")
print(search_res$ids)

# If it returns an ID, fetch summary
if (length(search_res$ids) > 0) {
  sum_res <- entrez_summary(db = "protein", id = search_res$ids[1])
  print(sum_res$taxid)
  print(sum_res$organism)
} else {
  message("Accession not found in protein DB")
}

#### write accessions to file to get taxonomy from dump ####
accessions_vec <- unlist(unique_ncbi_accessions)
writeLines(accessions_vec, "all_ncbi_accessions.txt")

accessions_vec <- unlist(unique_jgi_accessions)
writeLines(accessions_vec, "all_jgi_accessions.txt")
