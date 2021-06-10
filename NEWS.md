# MSnID 1.25


## Version 1.25.1

- Peptide mapping to protein sequence
- Visualization of peptides within protein

## Version 1.18.1

- Improved protein inference. Now it has two options - either unique peptide or
   with razor.

## Version 1.17.1

- Broadened the scope of the tool to ImmunoOncology - that is finding rare
  events that reflected in protein sequence.
   
## Version 1.15.1

- mzR now added scan number(s) into the table representation of
  mzIdentML object. As a results it caused an error in my unit test
   checking if the file reads properly. Fixed this check with updated hash.

## Version 1.11.1

- Switch from mclapply to parLapply in filter optimizaition.
  This fixes a bug 
  "Assertion failure at kmp_runtime.cpp(6480): __kmp_thread_pool == __null."
  Also allows to run parallel on Window OS

## Version 1.7.3

- Logical error fix in PSM FDR calculation. Prior it ignored redundancy
  in peptide to protein mapping. Now protein accession are discarded
  and unique rows considered before compuring the FDR.
- Added unit test to cover new PSM FDR calculation.

## Version 1.7.2

- added high performance parser based on mzR::openIDFile

## Version 1.7.1

- added infer_parsimonious_accessions method
- version leap from 1.3.1 to 1.7.1 to catch up with Bioconductor
- fixed conflict between reshape2 and data.table

## Version 1.3.1

- S4 methods "accessions" and "proteins" outsourced to ProtGenerics

## Version 1.1.6

- Added explicit mc.cores argument for "Grid" type of optimization
- Fix of trimming to basename on POSIX

## Version 1.1.5

- Always trimming spectrumFile and databaseFile to basename
  after reading data from mzIdentML files.

## Version 1.1.4

- Catching up with changes in mzID package version 1.5.3

## Version 1.1.3

- Updated data, documentation and unit test to reflect the
  latest changes in mzID

## Version 1.1.2

- id_quality and evaluate_filter may have multiple values for
  the level argument. By default it is all three:
  "PSM", "peptide" and "accession". The return value for
  both methods is switched to matrices.
- added Laurent Gatto as a contributor

## Version 1.1.1

- fix mishandling PSM
- use mzR::psms generic

## Version 1.1.0

- Bumping versions after creating 3.0 release branch.

## Version 0.99.3

- Rcpp back into "Depends"

## Version 0.99.2

- Updated Rcpp dependency
- Updated contact info to match the bioc-devel mailing list

## Version 0.99.1

- Updated dependency on mzID version

## Version 0.99.0

- Minor bug fixes
- Added unit tests for filter optimization and assessment of
  non-tryptic and missed cleavages
- Version bump for Bioconductor submission

## Version 0.1.0

- First release version
