library("MSnID")
library("dplyr")
library("Biostrings")

test_compute_accession_coverage <- function() {
  
  path_to_FASTA <- system.file("extdata", "for_phospho.fasta.gz", 
                               package = "MSnID")
  fasta <- readAAStringSet(path_to_FASTA)
  names(fasta) <- gsub(".*\\|(.*)\\|.*", "\\1", names(fasta))
  
  suppressMessages(msnid <- MSnID())
  psms(msnid) <- data.frame(
    accession = c(rep("B3KS81", 4), rep("O15075", 3)),
    peptide = c("MSSPKRSSKPSMSLAPS", "LAPSGSSMPTA", "TADPKPPASLKS", 
                "DNPSPSSSRK", "MSFGRDMELEHFDER", "FGRDME",
                "VKTTSASRAVSSLATAKGSPSEVR"),
    isDecoy = NA, calculatedMassToCharge = NA, 
    experimentalMassToCharge = NA, chargeState = NA, 
    spectrumFile = NA, spectrumID = NA)
  
  # Actual results
  out <- compute_accession_coverage(msnid, fasta, 
                                    accession_col = "accession", 
                                    pepSeq_col = "peptide") %>% 
    psms() %>% 
    select(accession, percentAACoverage) %>% 
    distinct()
  
  # Expected results
  ref <- data.frame(accession = c("B3KS81", "O15075"),
                    percentAACoverage = 100 * c(44/715, 39/740))
  
  checkIdentical(out, ref)
}


