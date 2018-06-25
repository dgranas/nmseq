# Nm-seq

Counts 5' ends of R2 reads from an Nm-seq experiment

## Prerequisites  
[samtools](http://www.htslib.org)  
python modules: pandas, (pysam is optional)

## Running on example data  
Click on 'Clone or Download' and select 'Download ZIP'  
In a terminal navigate to this directory `nmseq-master` and run:  
```python countEnds.py nm bam_sample.txt -r mouse_ref_seqs.fa -s shortnames.txt```  
Alternatively if the python pysam module is installed run:  
```python countEnds.py nm bam_sample.txt -r mouse_ref_seqs.fa -s shortnames.txt -p```  
If it works this should generate 4 .counts files 

## Running on new data  

### Create a bam_samples csv file  
Each line contains [bam filename], [sample name]  
The line order of the bam files determines the column order in the counts file    
This can be created manually or by running createBamSample.py

### Running countEnds.py  
python countEnds.py [experiment_name] [bam_samples file]
-r [reference sequences fasta]
-s [shortnames csv file]
-p

### Required arguments

* [experiment_name]  
name used in output files, for example nm2_rRNA or nm3_chrM

* [bam_samples file]  
csv file connecting bam filename to the sample name

### Optional arguments

* -r [reference sequences fasta file]  
fasta file with reference sequences  
headers must match bam refs (or the shortened names if using -s option)  
you would use this if you are analyzing rRNA alignments for example

* -s [shortnames csv file]  
csv file containing: [reference name], [shortened name]  
this is to convert the bam reference names to shorter names if desired  
for example 'gi|262231778|ref|NR_030686.1|' could be converted to '5S'

* -p
uses pysam module to directly read bam file instead of converting to sam  
use this command for slightly faster running time

