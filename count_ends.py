#!/usr/bin/env python3
'''
Given a file where each line has <bam filename>, <sample name>
For each bamfile:
Convert to sam file with R2 reads with SAM flag 131
(read paired, read mapped in proper pair, second in pair)
???Optional select only perfectly matched alignments???
Counts the 5' positions of the R2 reverse reads
'''

import argparse

def read_bamfiles_samples(bam_sample_files):
    '''
    Reads in file where each line has <bam file>, <sample name>
    Returns list of [(bam file, sample name)]
    '''
    # read in lines, avoiding any blank lines
    with open(bam_sample_files) as FH:
        lines = [line.rstrip() for line in FH.readlines() if line.rstrip()]

    # parse out the bam file and sample name which should be comma separated
    bamfiles_samples = []
    for line in lines:
        if ',' not in line:
            raise SystemExit('{0} does not contain a comma in line: {1}'.format(
                bam_sample_files, line))
        bamfile, sample = line.split(',')
        bamfiles_samples.append((bamfile, sample))
    return bamfiles_samples

def count_ends(bam_file, sample_name):
    '''
    Given a bam_file and sample_name
    Counts end
    '''

def main():
    parser = argparse.ArgumentParser(description= 'Count ends from bam alignments, requires samtools')
    parser.add_argument('bamfiles_samples', type=str, 
        help='file where each line: bam filename, sample name')
    parser.add_argument('-r', '--ref', type=str, 
        help='reference sequences in fasta format with headers that match bam references', 
        metavar='reference.fasta')
    args = parser.parse_args()

    bamfiles_samples = read_bamfiles_samples(args.bamfiles_samples)

    for bamfile, sample in bamfiles_samples:
        count_ends(bam_file, sample_name)

if __name__ == '__main__':
    main()