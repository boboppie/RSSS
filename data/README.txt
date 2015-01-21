This directory contains a copy of Ensembl DNA sequences, annotations, etc used to simulate the reads, and sample the transcript sets. They can be re-created by running:

Human:

wget -nH --cut-dirs=1 -r -l 1 -P ensembl_human/ ftp://ftp.ensembl.org/pub/release-66/fasta/homo_sapiens/dna/ &
wget -nH --cut-dirs=1 -r -l 1 -P ensembl_human/ ftp://ftp.ensembl.org/pub/release-66/fasta/homo_sapiens/cdna/ &
wget -nH --cut-dirs=1 -r -l 1 -P ensembl_human/ ftp://ftp.ensembl.org/pub/release-66/gtf/homo_sapiens/

Mouse:

wget -nH --cut-dirs=1 -r -l 1 -P ensembl_mouse/ ftp://ftp.ensembl.org/pub/release-66/fasta/mus_musculus/dna/ &
wget -nH --cut-dirs=1 -r -l 1 -P ensembl_mouse/ ftp://ftp.ensembl.org/pub/release-66/fasta/mus_musculus/cdna/ &
wget -nH --cut-dirs=1 -r -l 1 -P ensembl_mouse/ ftp://ftp.ensembl.org/pub/release-66/gtf/mus_musculus/

Worm:

wget -nH --cut-dirs=1 -r -l 1 -P ensembl_worm/ ftp://ftp.ensembl.org/pub/release-66/fasta/caenorhabditis_elegans/dna/ &
wget -nH --cut-dirs=1 -r -l 1 -P ensembl_worm/ ftp://ftp.ensembl.org/pub/release-66/fasta/caenorhabditis_elegans/cdna/ &
wget -nH --cut-dirs=1 -r -l 1 -P ensembl_worm/ ftp://ftp.ensembl.org/pub/release-66/gtf/caenorhabditis_elegans/

Yeast:

wget -nH --cut-dirs=1 -r -l 1 -P ensembl_yeast/ ftp://ftp.ensembl.org/pub/release-66/fasta/saccharomyces_cerevisiae/dna/ &
wget -nH --cut-dirs=1 -r -l 1 -P ensembl_yeast/ ftp://ftp.ensembl.org/pub/release-66/fasta/saccharomyces_cerevisiae/cdna/ &
wget -nH --cut-dirs=1 -r -l 1 -P ensembl_yeast/ ftp://ftp.ensembl.org/pub/release-66/gtf/saccharomyces_cerevisiae/