
map_peptides_to_fasta_file <- function(peptides, fasta, 
                                       numCores = detectCores() - 1) 
{

  fasta.df <- data.frame(accession = names(fasta),
                         sequence = fasta,
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

utils::globalVariables(c("is_match"))

