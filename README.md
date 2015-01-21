# A comparative study of RNA-seq analysis strategies
## user manual and guide

--------

## Overview

Three principal approaches have been proposed for inferring the set of transcripts expressed in RNA samples using RNA-seq. The simplest approach uses curated annotations, which assumes the transcripts in a sample are a subset of the transcripts listed in a curated database. A more ambitious method involves aligning reads to a reference genome and using the alignments to infer the transcript structures, possibly with the aid of a curated transcript database. The most challenging approach is to assemble reads into putative transcripts *de novo* without the aid of reference data. We have systematically assessed the properties of these three approaches through a simulation study. We have found that sensitivity of computational transcript set estimation is severely limited. Computational approaches (both genome-guided and *de novo* assembly) produce a large number of artefacts which are assigned large expression estimates and absorb a substantial proportion of the signal when performing expression analysis. The approach using curated annotations shows good expression correlation even when the annotations are incomplete. Furthermore, any incorrect transcripts present in a curated set do not absorb much signal, so it is preferable to have a curation set with high sensitivity than high precision. Software to simulate transcript sets, expression values and sequence reads under a wider range of parameter values and to compare sensitivity, precision and signal-to-noise ratios of different methods is freely available in this repository and can be expanded by interested parties to include methods other than the exemplars.

## Citing

A manuscript has been submitted to *Briefings in Bioinformatics* for peer review.

## Obtaining

To download the source code, please use git to download the most recent development
tree.  Currently, the tree is hosted on github, and can be obtained via:

    git clone git://github.com/boboppie/RSSS.git

## Usage

1. Download data (see data/README.txt for the details)
2. Install all the pre-required software
3. Set the **RSSS_DATA_DIR** environment variable to the path containing the Ensembl data
4. Execute `run_pipeline.sh` with the appropriate options: e.g. `run_pipeline.sh -c` to compute and visualise transcriptome reconstruction accuracy for varying coverage values, or `run_pipeline.sh -s` to compare sensitivity and precision values

## Contributors

RSSS is made by:

- Fengyuan Hu 
- Jürgen Jänes 
- Ernest Turro 
- Alexandra Lewin 

## Support

### email

Please report any issues or questions by creating a ticket, or by email to 
<fh293@cam.ac.uk>.
