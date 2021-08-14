parse_FASTA_names <- function(FASTA_names) {
  
  # Regex is limited to 9 groups, so gene, protein existence, 
  # and sequence version are grouped together and then split
  pttrn <- paste0("([a-z]{2})\\|([^-]+)\\-?([0-9]+)?\\|(.*_[A-Z]+)",
                  "\\s(.+?)\\sOS=(.*)\\sOX=(\\d+)\\s(.*)$")
  
  x <- data.frame(
    database = sub(pttrn, "\\1", FASTA_names),
    uniprot_acc = sub(pttrn, "\\2", FASTA_names),
    isoform = sub(pttrn, "\\3", FASTA_names),
    entry_name = sub(pttrn, "\\4", FASTA_names),
    description = sub(pttrn, "\\5", FASTA_names),
    organism = sub(pttrn, "\\6", FASTA_names),
    organism_id = sub(pttrn, "\\7", FASTA_names),
    gene = sub("GN=", "", 
               str_extract(sub(pttrn, "\\8", 
                               FASTA_names), "GN=[^\\s]+")),
    protein_existence = sub("PE=", "", 
                            str_extract(sub(pttrn, "\\8", 
                                            FASTA_names), "PE=\\d+")),
    sequence_version = sub("SV=", "", 
                           str_extract(sub(pttrn, "\\8", 
                                           FASTA_names), "SV=\\d+"))
  )
  
  return(x)
}

