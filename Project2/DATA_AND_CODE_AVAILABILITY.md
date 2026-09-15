# Data and code availability

## Data availability

This study uses individual-level HELIOS data and GWAS summary statistics from
external studies and consortia. These data are not redistributed through
GitHub. Access is governed by the relevant ethics approvals, data-use
agreements and provider application procedures described in the manuscript.

The repository contains no individual-level participant data. Non-sensitive
material supplied with the code includes analysis configurations, selected
derived matrices, run manifests, checksums, synthetic examples, archived
functional-enrichment exports and compact processed results. These files are
provided to identify the analyses reported in the paper and to support
verification of the computational workflow.

## Code availability

Custom R, Python, Bash and AWK code used for quality control, LDSC/GenomicSEM,
LAVA, PLEIO preparation and testing, directional locus classification,
SMR/HEIDI reanalysis, FUMA/MAGMA processing and g:Profiler consolidation is
available at:

<https://github.com/TimvanderEs/PhD/tree/main/Project2>

The repository includes workflow-specific instructions and automated tests
using synthetic data. Analyses that depend on controlled or third-party inputs
can be rerun after those inputs have been obtained from their original
providers. Third-party programs are not redistributed and should be installed
from their cited sources.

The version used for publication will be identified by a GitHub release and an
archival DOI. Until then, the relevant Git commit provides the immutable version
identifier.
