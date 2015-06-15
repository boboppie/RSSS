# A comparative study of RNA-seq analysis strategies

This repository holds the pipeline for simulating RNA-seq analyses published in:

> Jänes J, Hu F, Lewin A, and Turro E. [A Comparative Study of RNA-Seq Analysis Strategies.](http://dx.doi.org/10.1093/bib/bbv007) _Briefings in Bioinformatics_, March 2015. doi:10.1093/bib/bbv007.

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
