\name{parse_FASTA_names}
\alias{parse_FASTA_names}

\title{
Create table from UniProt FASTA headers
}

\description{
Separate the elements of UniProt FASTA headers into columns of a \code{data.frame}.
}

\usage{parse_FASTA_names(path_to_FASTA, database = c("uniprot", "gencode"))}

\arguments{
  \item{path_to_FASTA}{path to UniProt FASTA file.}
  \item{database}{database from which the FASTA file was derived. Currently accepts "uniprot" (default) or "gencode". Not case-sensitive.}
}

\value{
\code{data.frame} with database-specific columns for each component of the FASTA header.
}

\seealso{
\href{https://www.uniprot.org/help/fasta-headers}{UniProt FASTA headers},
\href{https://www.gencodegenes.org/pages/data_format.html}{GENCODE Data Format}
}

\author{
Tyler Sagendorf, Bennett Drucker
}

\examples{
path_to_FASTA <- system.file("extdata/uniprot_rat_small.fasta.gz",
                             package = "MSnID")
x <- parse_FASTA_names(path_to_FASTA)
head(x)
}
