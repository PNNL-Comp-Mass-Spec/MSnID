#' Map peptides to FASTA file
#'
#' @description
#' Takes list of clean peptide sequences and FASTA file and returns a
#' data frame containing all matches.
#' 
#' @param peptides (character) list of clean peptide sequences
#' @param fasta (AAStringSet)
#' @return (data.frame) with two columns for peptide and accession
#' @importFrom parallel detectCores makeCluster clusterEvalQ clusterExport clusterSplit clusterApply stopCluster
#' @importFrom dplyr %>% mutate filter select
#' 
#' @example 
#' msnid <- MSnID(".")
#' mzids <- system.file("extdata","phospho.mzid.gz",package="MSnID")
#' msnid <- read_mzIDs(msnid, mzids)
#' 
#' peptides <- unique(msnid$pepSeq)
#' 
#' fst_path <- system.file("extdata","for_phospho.fasta.gz",package="MSnID")
#' library(Biostrings)
#' fst <- readAAStringSet(fst_path)
#' 
#' out <- map_peptides_to_fasta_file(peptides, fst, numCores = 2L)


#' @export
map_peptides_to_fasta_file <- function(peptides, fasta, numCores = detectCores() - 1) {
  
  fasta.df <- data.frame(accession = names(fst),
                         sequence = fst,
                         row.names = NULL)
  
  cl <- makeCluster(getOption("cl.cores", numCores))
  clusterEvalQ(cl, library(dplyr))
  clusterExport(cl, "peptides")
  
  fasta.df.split <- lapply(clusterSplit(cl, 1:nrow(fasta.df)),
                          function(idx) {
                            fasta.df[idx,, drop=FALSE]
                          }) 
  
  f <- function(peptides, fasta.df) {
    lapply(peptides, function(peptide) {
      fasta.df %>%
        mutate(is_match = grepl(peptide, sequence)) %>%
        filter(is_match) %>%
        mutate(peptide = peptide) %>%
        select(peptide, accession, -sequence, -is_match)
      })
  }
  
  out <- clusterApply(cl, fasta.df.split, function(fasta.df) {
    f(peptides, fasta.df)})
  stopCluster(cl)
  
  out <- unlist(out, recursive = FALSE)
  out <- do.call(rbind, out)
}
