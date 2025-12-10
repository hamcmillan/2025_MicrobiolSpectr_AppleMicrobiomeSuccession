library(dplyr)
library(tidyr)
library(purrr)

# Define taxonomic levels from specific to general
tax_levels <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom")

# Function to get best match from CAZy
get_best_taxonomic_match <- function(taxon_row, cazy_df) {
  for (level in tax_levels) {
    query_value <- taxon_row[[level]]
    if (!is.na(query_value) && query_value != "") {
      # Check if there's a match at this level
      match_rows <- cazy_df %>%
        filter(!!sym(level) == query_value)
      if (nrow(match_rows) > 0) {
        # Return the first match (or sum if multiple)
        return(match_rows %>% 
                 summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
                 mutate(Match_Level = level))
      }
    }
  }
  # If no match found at any level, return NAs
  return(tibble(GH = NA, GT = NA, Match_Level = NA))
}

# Apply the function row-wise
df_results <- df_tax %>%
  rowwise() %>%
  mutate(match = list(get_best_taxonomic_match(cur_data(), cazy))) %>%
  unnest(cols = c(match))
