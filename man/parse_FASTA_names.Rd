\name{parse_FASTA_names}
\alias{parse_FASTA_names}

\title{
Create table from UniProt FASTA headers
}

\description{
Separate the elements of UniProt FASTA headers into columns of a \code{data.frame}.
}

\usage{parse_FASTA_names(path_to_FASTA)}

\arguments{
  \item{path_to_FASTA}{path to UniProt FASTA file.}
}

\value{
\code{data.frame} with the following columns: feature, database, uniprot_acc, isoform, entry_name, description, organism, organism_id, gene, protein_existence, and sequence_version.
}

\seealso{
\href{https://www.uniprot.org/help/fasta-headers}{UniProt FASTA headers}
}

\author{
Tyler Sagendorf
}

\examples{
path_to_FASTA <- system.file("extdata/uniprot_rat_small.fasta.gz",
                             package = "MSnID")
x <- parse_FASTA_names(path_to_FASTA)
head(x)
}
