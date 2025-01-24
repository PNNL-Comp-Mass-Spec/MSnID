\name{map_peptides_to_fasta_file}
\alias{map_peptides_to_fasta_file}
\title{Map peptides to FASTA file}

\usage{
map_peptides_to_fasta_file(peptides, 
                           fasta, 
                           numCores = detectCores() - 1)
}

\arguments{
\item{peptides}{(character) vector of clean peptide sequences.}

\item{fasta}{(AAStringSet object)}

\item{numCores}{(numeric) number of copies of \R to run in parallel.}
}

\value{
(data.frame) with two columns for peptide and accession
}

\description{
Takes list of clean peptide sequences and FASTA file and returns a
data frame containing all matches.
}

\examples{
msnid <- MSnID(".")
mzids <- system.file("extdata", "phospho.mzid.gz", package="MSnID")
msnid <- read_mzIDs(msnid, mzids)

peptides <- peptides(msnid)

path_to_FASTA <- system.file("extdata", "for_phospho.fasta.gz", 
                             package = "MSnID")
fasta <- Biostrings::readAAStringSet(path_to_FASTA)

out <- map_peptides_to_fasta_file(peptides, fasta, numCores = 2L)
}

\author{
Michael Nestor
}
