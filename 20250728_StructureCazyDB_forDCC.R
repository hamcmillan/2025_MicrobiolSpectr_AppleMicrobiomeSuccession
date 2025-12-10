install.packages(c("xml2", "data.table", "future.apply", "rentrez"))
library(xml2)
library(data.table)
library(future.apply)
library(rentrez)

# Set your registered NCBI API key here (increases limit to 10 req/sec)
rentrez::set_entrez_key("your_ncbi_api_key_here")
options(rentrez.delay = 0.12)  # 0.12s between requests = ~8.3 req/sec (safe)

get_taxonomy_from_accession_full <- function(accession) {
  tryCatch({
    search_res <- entrez_search(db = "protein", term = paste0(accession, "[ACCN]"))
    if (length(search_res$ids) == 0) {
      warning(paste("No protein ID found for", accession))
      return(NULL)
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
        return(NULL)
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
    
    data.frame(
      accession = accession,
      taxid = taxid,
      scientific_name = organism,
      t(major_lineage),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }, error = function(e) {
    warning(paste("Failed for accession:", accession, "-", e$message))
    return(NULL)
  })
}

# Load CAZy data
cazy <- read.table("cazy_data_20250707.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                   quote = "", fill = TRUE, comment.char = "",
                   col.names = c("Family", "Domain", "Organism", "Accession", "Source"))

# Filter only NCBI entries
cazy_ncbi <- cazy %>% filter(Source == "ncbi")

# Get unique accessions
unique_ncbi_accessions <- unique(cazy_ncbi$Accession)

# Your full list of accessions
accessions <- unique_ncbi_accessions  # length 3.8 million

# Split into chunks of 10,000
chunk_size <- 10000
chunks <- split(accessions, ceiling(seq_along(accessions) / chunk_size))

# Create a directory to save results
dir.create("taxonomy_batches", showWarnings = FALSE)

# Set up parallel plan
future::plan("multisession", workers = 10)  # Adjust based on DCC resource allocation

# Loop over chunks and process/save each in parallel
for (i in seq_along(chunks)) {
  chunk <- chunks[[i]]
  message("Processing chunk ", i, " of ", length(chunks))
  
  result <- future_lapply(chunk, get_taxonomy_from_accession_full)
  
  result_df <- rbindlist(result, fill = TRUE)
  fwrite(result_df, sprintf("taxonomy_batches/taxonomy_chunk_%04d.tsv", i), sep = "\t")
  
  # Optional: write log of failed accessions
  failed <- chunk[which(sapply(result, is.null))]
  writeLines(failed, sprintf("taxonomy_batches/failed_chunk_%04d.txt", i))
}

# Load all successfully completed chunks
file_list <- list.files("taxonomy_batches", pattern = "taxonomy_chunk_.*\\.tsv$", full.names = TRUE)
final_df <- rbindlist(lapply(file_list, fread), fill = TRUE)
fwrite(final_df, "cazy_data_20250707_taxonomy_full_results.tsv", sep = "\t")
