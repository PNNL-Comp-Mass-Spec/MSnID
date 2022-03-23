

parse_FASTA_names <- function(path_to_FASTA) {
  
  # Get FASTA headers
  fasta_names <- names(readAAStringSet(path_to_FASTA))
  
  # Pattern for first 6 columns
  pttrn <- "(.*)\\|(.*)\\|([^ ]+) (.*) OS=(.*) OX=(\\d+).*"
  
  x <- data.frame(feature = sub(" .*", "", fasta_names),
                  database = sub(pttrn, "\\1", fasta_names),
                  uniprot_acc = sub(pttrn, "\\2", fasta_names),
                  entry_name = sub(pttrn, "\\3", fasta_names),
                  description = sub(pttrn, "\\4", fasta_names),
                  organism = sub(pttrn, "\\5", fasta_names),
                  organism_id = sub(pttrn, "\\6", fasta_names),
                  # These last 3 may not be present, so we use separate regex
                  gene = sub(".* GN=([^ ]+).*", "\\1", fasta_names),
                  protein_existence = sub(".* PE=(\\d+).*", "\\1", fasta_names),
                  sequence_version = sub(".* SV=(\\d+)$", "\\1", fasta_names)
  ) %>% 
    mutate(isoform = str_extract(uniprot_acc, "(?<=-)\\d+"),
           uniprot_acc = sub("-.*", "", uniprot_acc),
           # By default, a failed match uses the original string. 
           # Replace with NA
           across(.cols = everything(),
                  .fns = ~ ifelse(.x == fasta_names, NA, .x))) %>% 
    # Remove entries that do not follow the format (no database entry)
    dplyr::filter(!is.na(database)) %>% 
    # Convert these columns from character to numeric
    mutate(across(.cols = c(isoform, protein_existence, 
                            sequence_version, organism_id),
                  .fns = as.numeric)) %>% 
    # Reorder columns
    dplyr::select(feature, database, uniprot_acc, isoform, everything())
  
  return(x)
}


utils::globalVariables(c("uniprot_acc", "database", "isoform", "feature",
                         "protein_existence", "sequence_version", 
                         "organism_id"))

