\name{infer_parsimonious_accessions}
\alias{infer_parsimonious_accessions}
\alias{infer_parsimonious_accessions,MSnID-method}

\title{Eliminates Redundancy in Peptide-to-Protein Mapping}

\description{
    Infer the parsimonious set of accessions (e.g. proteins) 
    that explains all of the peptide sequences. The algorithm is a 
    simple loop that looks for the accession explaining the most peptides,
    records the peptide-to-accession mapping for this accession,
    removes those peptides, and then looks for the next best accession. The
    loop continues until no peptides are left. This is known as the greedy set 
    cover algorithm.
}

\usage{
    infer_parsimonious_accessions(object, 
                                  unique_only=FALSE, 
                                  prior=character(0),
                                  refine_prior=FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{
        An instance of class "MSnID".
    }
    \item{unique_only}{
        If TRUE, peptides mapping to multiple accessions are dropped and
        only unique are retained. If \code{FALSE}, then shared peptides assigned
        according to Occam's razor rule. That is a shared peptide is 
        assigned to a protein with larger number of unique peptides.
        If the number of unique peptides is the same, then to the first 
        accession. Default is \code{FALSE}.
    }
    \item{prior}{(character)
        character vector with prior justified proteins/accessions.
        If \code{unique_only == TRUE}, then \code{prior} argument is
        ignored. Essentially, evidence by presence of unique peptide supercedes
        any prior. Default is \code{character(0)}: that is, no prior.
    }
    \item{refine_prior}{(logical)
        if \code{FALSE} (default), peptides are allowed to match multiple
        proteins in the prior. That is, the greedy set cover algorithm is only
        applied to the set of proteins not in the prior. If \code{TRUE}, the
        algorithm is applied to the prior and non-prior sets separately before
        combining.
    }
}

\details{
    Although the algorithm is rather simple, it is \emph{THE} algorithm used for
    inferring maximal matching in bipartite graphs and is used in the 
    IDPicker software (see References).
}

\value{
    Returns an instance of "MSnID" with minimal set of proteins
    necessary to explain all the peptide sequences.
}

\author{
    Vladislav A Petyuk \email{vladislav.petyuk@pnnl.gov},
    Colin Brislawn,
    Tyler Sagendorf
}

\seealso{
    \code{\link{MSnID}}
}

\references{
Zhang, B., Chambers, M. C., & Tabb, D. L. (2007). Proteomic parsimony through bipartite graph analysis improves accuracy and transparency. Journal of Proteome Research, 6(9), 3549-3557. PubMed. \url{https://doi.org/10.1021/pr070230d}
}


\examples{
data(c_elegans)

# explicitely adding parameters that will be used for data filtering
msnidObj$msmsScore <- -log10(msnidObj$`MS-GF:SpecEValue`)
msnidObj$absParentMassErrorPPM <- abs(mass_measurement_error(msnidObj))

# quick-and-dirty filter. The filter is too strong for the sake of saving time
# at the minimal set of proteins inference step.
msnidObj <- apply_filter(msnidObj, 'msmsScore > 12 & absParentMassErrorPPM < 2')

show(msnidObj)
msnidObj2 <- infer_parsimonious_accessions(msnidObj)
show(msnidObj2)

}


