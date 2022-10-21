

parse_FASTA_names <- function(path_to_FASTA,
                              database = c("uniprot", "gencode"))
{
  database <- match.arg(arg = tolower(database),
                        choices = c("uniprot", "gencode"))
  
  # Get FASTA headers
  fasta_names <- names(readAAStringSet(path_to_FASTA))
  
  if (database == "uniprot") {
    # Pattern for first 6 columns
    pttrn <- "(.*)\\|(.*)\\|([^ ]+) (.*) OS=(.*) OX=(\\d+).*"
    
    out <- data.frame(
      feature = sub(" .*", "", fasta_names),
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
      mutate(
        isoform = str_extract(uniprot_acc, "(?<=-)\\d+"),
        uniprot_acc = sub("-.*", "", uniprot_acc),
        # By default, a failed match uses the original string.
        # Replace with NA
        across(
          .cols = everything(),
          .fns = ~ ifelse(.x == fasta_names, NA, .x)
        )
      ) %>%
      # Remove entries that do not follow the format (no database entry)
      dplyr::filter(!is.na(database)) %>%
      # Convert these columns from character to numeric
      mutate(across(
        .cols = c(
          isoform, protein_existence,
          sequence_version, organism_id
        ),
        .fns = as.numeric
      )) %>%
      # Reorder columns
      dplyr::select(feature, database, uniprot_acc, isoform, everything())
    
  } else if (database == "gencode") {
    pttrn <- paste(rep("([^\\|]+)", times = 8), collapse = "\\|")
    pttrn <- paste0(pttrn, "(.*$)")
    
    out <- data.table(value = fasta_names)
    out[, `:=`(
      protein_id = sub(pttrn, "\\1", value),
      transcript_id = sub(pttrn, "\\2", value),
      gene_id = sub(pttrn, "\\3", value),
      havana_gene = sub(pttrn, "\\4", value),
      havana_transcript = sub(pttrn, "\\5", value),
      transcript = sub(pttrn, "\\6", value),
      gene = sub(pttrn, "\\7", value),
      prot_length = as.numeric(sub(pttrn, "\\8", value)),
      extra_info = sub(pttrn, "\\9", value),
      value = NULL
    )]
    # Clean up columns
    out[out == "" | out == "-"] <- NA
    # Remove columns with all missing values
    out <- out[, which(unlist(
      lapply(out, function(x) !all(is.na(x)))
    )), with = FALSE]
    out <- as.data.frame(out)
  }
  
  return(out)
}


utils::globalVariables(
  c(
    "uniprot_acc",
    "database",
    "isoform",
    "feature",
    "protein_existence",
    "sequence_version",
    "organism_id",
    "value"
  )
)

