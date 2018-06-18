# Nm-seq

Counts 5' ends of R2 reads from an Nm-seq experiment

## Getting Started

### Prerequisites

[samtools](http://www.htslib.org)

python modules: pandas, argparse

### Create a bam_samples csv file  
Each line contains [bam filename], [sample name]  
This can be created manually or by running CreateBamSample.py

### Running CountEnds.py  
python CountEnds.py [experiment_name] [bam_samples file]  
-r [reference sequences fasta file]  
-s [shortnames csv file]

[experiment_name]  
name used in output files, for example nm2_rRNA or nm3_chrM

[bam_samples file]  
csv file connecting bam filename to the sample name

-r [reference sequences fasta file]  
fasta file with reference sequences  
headers must match bam refs  
you would use this if you are analyzing rRNA alignments for example

-s [shortnames csv file]  
csv file containing: [reference name], [shortened name]  
this is to convert the bam reference names to shorter names if desired  
for example 'gi|262231778|ref|NR_030686.1|' could be converted to '5S'

