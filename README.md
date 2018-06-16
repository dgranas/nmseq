# Nm-seq

Count ends from an Nm-seq experiment

## Getting Started

### Prerequisites

samtools needs to be installed

### Running the program

python count_ends.py [bamfiles_samples] [output]

[bamfiles_samples] is a csv text file containing:
bame filename, sample name

[output] is the name for the output, for example nm2_rRNA

Optional arguments:

-r, --ref: fasta file with reference sequences, headers must match bam refs

-mm, --mismatch: allow mismatches in bam alignments (default is no mismatches)

-f, --samflag: set the [SAMFLAG](https://broadinstitute.github.io/picard/explain-flags.html)
  default is 131 for read paired, read mapped in proper pair, second in pair

--header: keep header in sam file (default is no header)

--keepsam: keep sam files (default is they are deleted)

--refnames: csv file containing: [reference name], [shortened name]
  for example: gi|12044..., 28S