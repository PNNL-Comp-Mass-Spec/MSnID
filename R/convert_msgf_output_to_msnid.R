
utils::globalVariables(c("accession", "Protein",
                         "calculatedMassToCharge", "MH", "Charge",
                         "chargeState", "experimentalMassToCharge",
                         "PrecursorMZ", "isDecoy",
                         "spectrumFile", "Dataset",
                         "spectrumID", "Scan",
                         "peptide", "Peptide",
                         "pepSeq"))

convert_msgf_output_to_msnid <- function(x) {
  suppressMessages(msnid <- MSnID("."))
  x <- x %>% mutate(accession = Protein,
                    calculatedMassToCharge = (MH + (Charge-1)*.PROTON_MASS)/Charge,
                    chargeState = Charge,
                    experimentalMassToCharge = PrecursorMZ,
                    isDecoy = grepl("^XXX", Protein),
                    spectrumFile = Dataset,
                    spectrumID = Scan) %>%
    rename(peptide = Peptide)
  # clean peptide sequences
  x <- mutate(x, pepSeq = .get_clean_peptide_sequence(peptide))
  
  psms(msnid) <- x
  
  return(msnid)
}