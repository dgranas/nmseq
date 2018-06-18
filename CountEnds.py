#!/usr/bin/env python3
'''
Requires: samtools

python CountEnds.py [bamfiles_samples] [experiment_name]

[bamfiles_samples]: csv file where each line has <bam filename>, <sample name>
[experiment_name]: name used in output files, for example nm2_rRNA or nm3_chrM 

optional arguments:
-r, --ref: fasta file with reference sequences, headers must match bam refs
--refnames: csv file containing: [reference name], [shortened name]

For each bamfile:
Convert to sam file with R2 reads with SAM flag 131
(read paired, read mapped in proper pair, second in pair)
???Optional select only perfectly matched alignments???
Counts the 5' positions of the R2 reverse reads
'''

import argparse
import os
import collections
import pandas as pd

def fasta_to_dict(fasta_file):
    '''
    Given a fasta file
    Returns a dict with {header: sequence}
    '''
    header_to_seq = {}
    with open(fasta_file) as FH:
        for line in FH:
            line = line.rstrip()
            if line.startswith('>'):
                header = line[1:]
                header_to_seq[header] = ''
            else:
                header_to_seq[header] += line.upper()
    return header_to_seq

def bam_to_sam(bamfile, sam_flag=131, mismatch=False, keep_header=False):
    '''
    Given a bam file
    Converts it to a sam file 
    Applies -f sam flag 
    For example to get only R2 reads you can use -f 131
    131 = read paired, read mapped in proper pair, second in pair
    If mismatch is False it applies an awk command to select for NM:i:0
    If keep_header is True is uses -h to keep the header
    Returns the sam file
    '''
    # change <file>.bam to <file>.sam
    file_prefix = os.path.join(*os.path.splitext(bamfile)[:-1])
    samfile = '{}.f{}.sam'.format(file_prefix, sam_flag)

    command = 'samtools view {} -f {}'.format(bamfile, sam_flag)

    if keep_header:
        command += ' -h'
    if not mismatch:
        command += " | awk '$0 ~\"NM:i:0\"'"

    command += ' > {}'.format(samfile)

    os.system(command)
    return samfile   

def count_ends(samfile, ref_to_seq, ref_to_shortname):
    '''
    Given a sam file, reference sequences, and shortnames
    Count the 5' ends of reads
    '''
    counts = {}

    with open(samfile) as FH:
        for line in FH:
            line = line.rstrip()

            l = line.split()
            flag = int(l[1])
            ref_name = l[2]
            pos = int(l[3]) # left-most 1-based position
            cigar = l[5]
            seq = l[9]

            # assuming no mismatches allowed, all cigars should be [int]M
            ref_length = int(cigar[:-1])

            if flag & 16: # read is reverse strand, so 5' is the right-most position
                pos = pos + ref_length - 1
                # -1 since if starts at pos 1 and is 3bp long, it will go to pos 1+3-1 = 3
                read_base = seq[-1]
            else: # read is forward strand, so 5' is the left-most base
                read_base = seq[0] 

            ref_name = ref_to_shortname.get(ref_name, ref_name)
            if ref_to_seq:
                # check if the ref_name doesn't match any header in the ref fasta file
                if ref_name not in ref_to_seq.keys():
                    raise SystemExit('{} in sam file did not match fasta headers'.format(ref_name))

            if ref_name not in counts.keys():
                counts[ref_name] = collections.defaultdict(int)

            counts[ref_name][pos] += 1
    return counts

def generate_rows(seen_pos, ref_range=None):
    '''
    Given a set of seen positions
    If ref_range is None (genome case) it yields just the seen positions
    Else it yields all positions in from 0 to ref_range-1
    '''
    if not ref_range: # genome case
        for pos in sorted(seen_pos):
            yield pos
    else: # rRNA case
        for pos in range(ref_range):
            yield pos

def write_counts(counts, output, samples, ref_to_seq, ref_to_shortname):
    '''
    Writes end counts to file, each col is a sample, each row is a position
    If ref_to_seq is not None then it will:
        -output a column for bases from the ref seqs
        -output counts for each position even if 0 in all samples
        -create separate files for each reference sequence
    If ref_to-seq is None it will only output counts seen in any sample
    '''
    # get all the ref_names that had counts
    ref_names = set()
    for sample in samples:
        for ref_name in counts[sample].keys():
            ref_names.add(ref_name)

    # get all positions that were seen in at least 1 sample
    seen_pos = {}
    if not ref_to_seq: # genome case 
        for ref_name in ref_names:
            seen_pos[ref_name] = set()
            for sample in samples:
                seen_pos[ref_name].update(counts[sample][ref_name].keys())
        
        #outfile = f'{output}.counts'
        outfile = '{}.counts'.format(output)
        ref_range = None # what is ref_range??????
    
    file_opened = False

    for ref_name in sorted(ref_names):
        # generate the header for each ref_name
        header = ['ref_name', 'pos']
        if ref_to_seq:
            header.append('base')
            # for rrna case write to separate files for each ref_name
            #outfile = f'{output}_{ref_name}.counts'
            outfile = '{}_{}.counts'.format(output, ref_name)
            ref_range = len(ref_to_seq[ref_name])
        for sample in samples:
            header.append(sample)

        if os.path.exists(outfile) and not file_opened:
            raise SystemExit('file {} already exists'.format(outfile))

        with open(outfile, 'a') as f:
            f.write('\t'.join(header) + '\n')    
            for pos in generate_rows(seen_pos.get(ref_name, None), ref_range):
                row = [ref_name, pos]
                # add the corresponding base if reference seq was given
                if ref_to_seq:
                    row.append(ref_to_seq[ref_name][pos])
                for sample in samples:
                    row.append(counts[sample][ref_name][pos])
                f.write('\t'.join([str(i) for i in row]) + '\n')
            if not ref_to_seq:
                file_opened = True

def main():
    parser = argparse.ArgumentParser(
        description= 'Count ends from bam alignments, requires samtools')
    parser.add_argument('bamfiles_samples', type=str, 
        help='file where each line: bam filename, sample name')
    parser.add_argument('experiment_name', type=str,
        help='experiment name (for example nm2_rRNA, or nm3_chrM)')
    parser.add_argument('-r', '--refseqs', type=str, 
        help='reference sequences in fasta format with headers that match bam references', 
        metavar='reference.fasta')
    parser.add_argument('-s', '--shortnames', type=str,
        help='file containing <reference name>, <shortened name> (gi|12044..., 28S)')
    args = parser.parse_args()

    # if reference sequences provided, read them in and store as a dict {header: seq}
    if args.refseqs:
        args.refseqs = fasta_to_dict(args.refseqs)
        # make sure ref seqs start with hyphen to be 0-based to match sam positions
        for header, seq in args.refseqs.items():
            if seq[0].upper() in 'ACGT':
                args.refseqs[header] = '-' + seq

    # if shortened ref names are provided, get a dict for {ref: shortname}
    ref_to_shortname = {}
    if args.shortnames:
        shortnames = pd.read_csv(args.shortnames, names=['ref', 'shortname'], skipinitialspace=True)
        ref_to_shortname = dict(zip(shortnames.ref, shortnames.shortname))

    # read in bamfiles and create dict of {bamfile: sample}
    bamfiles_samples = pd.read_csv(args.bamfiles_samples, names=['bamfile', 'samples'],
        skipinitialspace=True)
    bamfile_to_sample = dict(zip(bamfiles_samples.bamfile, bamfiles_samples.samples))

    counts = {}
    for bamfile, sample in bamfile_to_sample.items():

        # convert bam file to sam file
        samfile = bam_to_sam(bamfile)

        # do the end counting
        counts[sample] = count_ends(samfile, args.refseqs, ref_to_shortname)

        # remove sam file
        if os.path.exists(samfile):
            os.remove(samfile)

    # write end counts to file
    write_counts(counts, args.experiment_name, bamfiles_samples.samples, args.refseqs,
                 ref_to_shortname)

if __name__ == '__main__':
    main()