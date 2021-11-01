\name{parse_FASTA_names}
\alias{parse_FASTA_names}

\title{
Create table from UniProt FASTA headers
}

\description{
Separate elements of UniProt FASTA names into columns of a data frame.
}

\usage{parse_FASTA_names(FASTA_names)}

\arguments{
  \item{FASTA_names}{UniProt FASTA headers}
}

\value{
data frame with the following columns: database, uniprot_acc, isoform, entry_name, description, organism, organism_id, gene, protein_existence, and sequence_version.
}

\note{
See \href{https://www.uniprot.org/help/fasta-headers}{https://www.uniprot.org/help/fasta-headers} for more info.
}

\examples{
\dontrun{
library(Biostrings)
# Need path_to_FASTA
fst <- readAAStringSet(path_to_FASTA)
FASTA_tbl <- parse_FASTA_names(names(fst))
}
}

