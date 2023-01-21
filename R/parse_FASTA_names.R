

parse_FASTA_names <- function(path_to_FASTA,
                              database = c("uniprot", "gencode"))
{
  database <- match.arg(arg = tolower(database),
                        choices = c("uniprot", "gencode"))
  
  # Get FASTA headers
  fasta_names <- names(readAAStringSet(path_to_FASTA))
  
  # Remove contaminants
  fasta_names <- fasta_names[!grepl("^contaminant", fasta_names, 
                                    ignore.case = TRUE)]
  
  if (database == "uniprot") {
    pttrn <- "(?<header>(?<feature>^(?<database>[^\\|]+)\\|(?<unique_id>(?<uniprot_acc>[^\\|-]+)-?(?<isoform>\\d*))\\|(?<entry_name>\\w+))\\s+(?<description>.*)\\s+OS=(?<organism>.*)\\s+OX=(?<organism_id>\\d+).*$)"
    
    m <- regexec(pattern = pttrn, text = fasta_names, perl = TRUE)
    out <- do.call(rbind, lapply(regmatches(fasta_names, m), `[`, -1L))
    out <- as.data.table(out)
    
    # Additional columns 
    # (Including positive lookbehind in pttrn removes isoform rows...)
    out[, `:=`(
      gene = sub(".*(?<=GN=)(\\S+).*", "\\1", header, perl = TRUE),
      protein_existence = sub(".*(?<=PE=)(\\d+).*", "\\1", header, perl = TRUE),
      sequence_version = sub(".*(?<=SV=)(\\d+)", "\\1", header, perl = TRUE),
      header = NULL
    )]
    
    # Convert these columns to numeric
    num_cols <- c("isoform", "protein_existence", 
                  "sequence_version", "organism_id")
    
  } else if (database == "gencode") {
    # Expect 8 sections separated by pipe
    pttrn <- paste(rep("([^\\|]+)", times = 8), collapse = "\\|")
    pttrn <- paste0(pttrn, ".*") # account for any additional information
    
    # Much faster than calling sub() separately for each column
    m <- regexec(pattern = pttrn, text = fasta_names)
    out <- do.call(rbind, lapply(regmatches(fasta_names, m), `[`, -1L))
    out <- as.data.table(out)
    colnames(out) <- c("protein_id", "transcript_id", "gene_id", "havana_gene", 
                       "havana_transcript", "transcript", "gene", 
                       "protein_length")
    
    # Convert these columns to numeric
    num_cols <- c("protein_length")
  }
  
  # Replace failed matches with NA
  for (j in seq_along(out)) {
    set(out, i = which(out[[j]] %in% c(fasta_names, "-", "")), j = j, 
        value = NA)
  }
  
  # Remove rows with all missing values
  out <- out[rowSums(!is.na(out)) != 0, ]
  
  # Convert certain columns to numeric
  out[, num_cols] <- out[, lapply(.SD, as.numeric), .SDcols = num_cols]
  
  # Remove columns with all missing values
  out <- out[, which(unlist(
    lapply(out, function(x) { !all(is.na(x)) })
  )), with = FALSE]
  
  setDF(out)
  
  return(out)
}


utils::globalVariables(
  c("header")
)

