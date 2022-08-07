\name{plot_protein_coverage}
\alias{plot_protein_coverage}
\alias{plot_protein_coverage,MSnID,character-method}
\title{Shows the peptide location(s) on the protein sequence}

\description{
Shows the peptide location(s) on the protein sequence.
}

\usage{
    plot_protein_coverage(object, 
                          accession, 
                          peptide_fill = "spectral.counts",
                          name_from = c("accession"),
                          border_size = NULL,
                          border_color = "white",
                          aa_step = NULL,
                          save_plot = FALSE,
                          ...)
}


\arguments{
    \item{object}{
        An instance of class MSnID.
    }
    \item{accession}{
        (string) accession name.
    }
    \item{peptide_fill}{
        (string) "spectral.counts", "sample.counts", or a column from 
        \code{psms(object)} used to fill the peptide bars.
    }
    \item{name_from}{
        (string) column sources for the plot title.
    }
    \item{border_size}{
        (numeric) controls the border around the peptide rectangles.
        In case of NULL, computes the value based on the protein length.
    }
    \item{border_color}{
        (character) color around peptide rectangles.
    }
    \item{aa_step}{
        (numeric) interval for the tickmarks along the protein sequence.
        In case of default (NULL) the tickmars are precomputed using 
        \code{breaks_pretty(30)}
    }
%    \item{...}{
%        Used for passing save_plot. If \code{TRUE}, then image is saved as a 
%        PNG file.
%    }
    \item{save_plot}{
        (boolean) if TRUE, the plot is saved as png file.
    }
}





\value{
    if save_plot = FALSE, then returns ggplot2 object
}

\author{
    Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov}
}


\examples{
m <- MSnID(".")
mzids <- system.file("extdata","phospho.mzid.gz",package="MSnID")
m <- read_mzIDs(m, mzids)
fst_path <- system.file("extdata","for_phospho.fasta.gz",package="MSnID")
library(Biostrings)
fst <- readAAStringSet(fst_path)
# ensure names are the same format as accessions(m)
names(fst) <- sub("(^[^ ]*) .*$", "\\1", names(fst))

# Creating sequences with repeats. This is just for the sake of 
# demonstration of the capability.
for(i in 2:4){
    fst[i] <- paste0(as.character(fst[i]),as.character(fst[i]))
}

# Appending reverse hits. Also, just for the demonstrating the capability.
mr <- mf <- apply_filter(m, "!isDecoy")
library(stringi)
mr$peptide <- stringi::stri_reverse(mr$peptide)
mr$accession <- paste0("XXX_", mr$accession)
mr$isDecoy <- TRUE
psms(m) <- rbind(psms(mr), psms(mf))

# the main call
m2 <- map_peptide_position(m, fst)
plot_protein_coverage(m2, "sp|Q9Y2W1|TR150_HUMAN", save_plot = FALSE)

# clean-up cache
unlink(".Rcache", recursive=TRUE)
}

