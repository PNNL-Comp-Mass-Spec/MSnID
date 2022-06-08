\name{compute_accession_coverage}
\alias{compute_accession_coverage}
\alias{compute_accession_coverage,MSnID-method}
\title{Compute accession coverage from identified peptides}

\description{
Given a list of accession and peptide pairs, computes the proportion of amino acid sequence covered by the listed peptides for each accession.
}

\usage{
    compute_accession_coverage(object,
                               fasta,
                               accession_col="accession",
                               pepSeq_col="pepSeq",
                               remove_nonmapping_peptides = TRUE)
}


\arguments{
    \item{object}{
        An instance of class MSnID.
    }
    \item{fasta}{
        (AAStringSet object) Protein sequences read from a FASTA file.
        Names must match the accessions of the MSnID object.
    }
    \item{accession_col}{
        (character) Name of the accession column.
    }
    \item{pepSeq_col}{
        (character) Name of the peptide sequence column.
    }
    \item{remove_nonmapping_peptides}{
        (boolean) Default it TRUE. If set to FALSE and non-mapped peptides are 
        present, there will be an error.
    }
}

\value{
An object of class MSnID with additional column `percentAACoverage` that contains values in 0-100.
}

\author{
Tyler Sagendorf, Michael Nestor, Vlad Petyuk
}


\examples{
# Create MSnID object
m <- MSnID(".")
mzids <- system.file("extdata", "phospho.mzid.gz", package = "MSnID")
m <- read_mzIDs(m, mzids)

# Get protein sequences
fst_path <- system.file("extdata", "for_phospho.fasta.gz", package = "MSnID")
fst <- Biostrings::readAAStringSet(fst_path)

# Remove decoys because they are not in this small FASTA file
m <- apply_filter(m, "!isDecoy")

# Compute accession coverage
m <- compute_accession_coverage(m, fst)

head(unique(subset(psms(m), select = c("accession", "percentAACoverage"))))

# Clean-up cache
unlink(".Rcache", recursive=TRUE)
}


