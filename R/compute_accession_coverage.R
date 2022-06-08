
.compute_accession_coverage <- function(object,
                                        fasta,
                                        accession_col,
                                        pepSeq_col,
                                        remove_nonmapping_peptides) 
{
  # Clean up on exit - optional, but frees up memory
  on.exit(rm(list = ls()))
  on.exit(gc(verbose = FALSE), add = TRUE)
  
  ## MS/MS data.table ----
  msms_cols <- c(accession_col, pepSeq_col)
  msms.dt <- unique(object@psms[, msms_cols, with = FALSE])
  
  ## Protein sequence data.table ----
  # Include reverse sequences if there are decoys
  decoy_acc <- unique(apply_filter(object, "isDecoy")[[accession_col]])
  if (length(decoy_acc) > 0 & !any(decoy_acc %in% names(fasta))) {
    fasta_rev <- reverse(fasta)
    names(fasta_rev) <- paste0("XXX_", names(fasta))
    fasta <- c(fasta, fasta_rev)
  }
  
  # check if fasta entry names are unique
  if(any(duplicated(names(fasta)))) {
    stop("FASTA entry names are not unique!\n")
  }
  
  # check if there is at least some agreement in IDs
  if (!all(object[[accession_col]] %in% names(fasta))) {
    stop("Some accession IDs not found in FASTA entries!\n")
  }
  
  # Reduce fasta to relevant accessions and create data.table
  fasta <- fasta[names(fasta) %in% msms.dt[[accession_col]]]
  fasta.dt <- data.table(width = width(fasta), 
                         seq = as.character(fasta),
                         accession_col = names(fasta))
  
  # Merge MSMS and FASTA data.tables
  msms.dt <- merge(msms.dt, fasta.dt, 
                   by.x = accession_col, by.y = "accession_col")
  
  ## Calculate overlap ----
  # Get starting positions of each peptide in corresponding
  # protein sequence by accession and calculate width of each peptide
  msms.dt[, start := as.numeric(stringr::str_locate(pattern = get(pepSeq_col),
                                                    string = seq)[, 1]),
          by = accession_col][, peptide_width := nchar(get(pepSeq_col))]
  
  # Safety feature. At this point we aren't sure how to handle peptides 
  # that are missing besides removing them.
  if(remove_nonmapping_peptides)
     msms.dt <- msms.dt[!is.na(start)]
  
  # Create IRanges object to calculate overlap
  irl <- IRanges(start = msms.dt[, start],
                 width = msms.dt[, peptide_width],
                 names = msms.dt[, pepSeq_col, with = FALSE])
  irl <- split(irl, msms.dt[, accession_col, with = FALSE]) # convert to list
  
  # Calculate AA overlap - this part takes the longest
  fullAACoverage <- sapply(irl, function(irl_i) {
    sum(IRanges::width(reduce(irl_i)))
  })
  # Convert named list to data.table
  fullAACoverage.dt <- data.table(fullAACoverage = fullAACoverage,
                                  accession_col = names(fullAACoverage))
  
  # Add fullAACoverage column to dt
  msms.dt <- merge(msms.dt, fullAACoverage.dt,
                   by.x = accession_col, by.y = "accession_col")
  # Calculate percentAACoverage and only keep required columns
  msms.dt <- unique(
    msms.dt[, percentAACoverage :=
              100 * fullAACoverage / width][
                , c(accession_col, "percentAACoverage"), with = FALSE]
  )
  
  # Add percentAACoverage column to psms(object)
  object@psms <- merge(object@psms, msms.dt, by = accession_col)
  return(object)
}

utils::globalVariables(c("start", "peptide_width", "percentAACoverage"))



