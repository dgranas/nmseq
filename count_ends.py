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
import os
import re
import collections
import pandas as pd

def get_cigar_len(newcig):
    '''
    Given a cigar string
    Returns the length of the alignment on the reference sequence
    '''
    cigar_len = 0

    while(len(newcig) > 0):
        ref_align = re.match(r'(\d+)[MNDP\=X](.*)', newcig)
        ref_unalign = re.match(r'(\d+)[SIH](.*)', newcig)

        if ref_align:
            cigar_len = cigar_len + int(ref_align.group(1))
            newcig = ref_align.group(2)
        elif ref_unalign:
            newcig = ref_unalign.group(2)
        else:
            raise SystemExit('Failed to parse CIGAR {}'.format(newcig))
            newcig = ''
    return cigar_len

def fasta_to_dict(fasta_file):
    '''
    Given a fasta file
    Returns a dict with {header: sequence}
    '''
    header_to_seq = {}
    header = None
    with open(fasta_file) as FH:
        for line in FH:
            line = line.rstrip()
            if line.startswith('>'):
                header = line[1:]
                header_to_seq[header] = ''
            else:
                header_to_seq[header] += line
    return header_to_seq

def read_bamfiles_samples(bamfiles_samples):
    '''
    Reads in file where each line has <bam file>, <sample name>
    Returns list of [(bam file, sample name)]
    '''
    # read in lines, avoiding any blank lines
    with open(bamfiles_samples) as FH:
        lines = [line.rstrip() for line in FH.readlines() if line.rstrip()]

    # parse out the bam file and sample name which should be comma separated
    bamfiles_samples = []
    for line in lines:
        if ',' not in line:
            raise SystemExit('{0} does not contain a comma in line: {1}'.format(
                bamfiels_samples, line))
        bamfile, sample = line.split(',')
        bamfiles_samples.append((bamfile.strip(), sample.strip()))
    return bamfiles_samples

def read_two_column_file

def bam_to_sam(bamfiles_samples, sam_flag, no_mismatch=True, keep_header=False):
    '''
    Given list of [(bamfile, sample)]
    Converts each bam file to a sam file 
    Applies -f sam flag 
    For example to get only R2 reads you can use -f 131
    If no_mismatch is True it applies an awk command to select for NM:i:0
    If keep_header is True is uses -h to keep the header
    Returns a list of [(sam_file, sample)]
    '''
    samfiles_samples = []
    for bamfile, sample in bamfiles_samples:
        # change <file>.bam to <file>.sam
        file_prefix = os.path.join(*os.path.splitext(bamfile)[:-1])
        samfile = '{}.f{}.sam'.format(file_prefix, sam_flag)

        command = 'samtools view {0} -f {1}'.format(bamfile, sam_flag)

        if keep_header:
            command += ' -h'

        if no_mismatch:
            command += " | awk '$0 ~\"NM:i:0\"'"

        command += ' > {}'.format(samfile)
        os.system(command)
        samfiles_samples.append((samfile, sample))
    return samfiles_samples

def count_ends(samfile, sample, ref_to_seq):
    '''
    Given a samfile and sample and reference sequences
    Counts the 5' ends of reads
    '''
    counts = {}

    with open(samfile) as FH:
        for line in FH:
            line = line.rstrip()

            if line.startswith('@'):
                continue

            l = line.split()
            flag = int(l[1])
            ref_name = l[2]
            pos = int(l[3]) # left-most 1-based position
            cigar = l[5]
            seq = l[9]
            
            if ref_to_seq:
                # check if the ref_name doesn't match any header in the ref fasta file
                if ref_name not in ref_to_seq.keys():
                    raise SystemExit('{} in sam file does not match any fasta headers'.format(
                        ref_name))
                # check if the ref base matches what the ref seq says
                if read_base != ref_to_seq[ref_name][pos]:
                    raise SystemExit('{} does not match {} at pos {} in {}'.format(
                        read_base, ref_to_seq[ref_name][pos], pos, ref_name))

            if ref_name not in counts.keys():
                counts[ref_name] = collections.defaultdict(int)

            if flag & 16: # read is reverse strand, so 5' is the right-most position
                pos = pos + get_cigar_len(cigar) - 1
                # -1 since if starts at pos 1 and is 3bp long, it will go to pos 1+3-1 = 3
                read_base = seq[-1] 
            else: # read is forward strand, so 5' is the left-most base
                read_base = seq[0]

            counts[ref_name][pos] += 1

    return counts

def generate_rows(seen_pos, ref_range=None):
    '''
    Given a counts dict for a particular ref_name
    if ref_range is None it will generate all the positions in the counts_ref
    else it will generate all positions in range(ref_range)
    '''
    if not ref_range: # genome case
        for pos in sorted(seen_pos):
            yield pos
    else: # rRNA case
        for pos in range(ref_range):
            yield pos

def write_counts(counts, output, samples, ref_to_seq):
    '''
    Writes end counts to file, each col is a sample, each row is a position
    If ref_to_seq is not None then it will:
        -output a column for bases from the ref seqs
        -output counts for each position even if 0 in all samples
        -create separate files for each reference sequence
    If ref_to-seq is None it will only output counts seen in any sample
    '''
    # get all the ref_names
    ref_names = set()
    for sample in samples:
        for ref_name in counts[sample].keys():
            ref_names.add(ref_name)

    # get all positions that were seen in at least 1 sample
    if not ref_to_seq: # genome case 
        seen_pos = {}
        for ref_name in ref_names:
            seen_pos[ref_name] = set()
            for sample in samples:
                seen_pos[ref_name].update(counts[sample][ref_name].keys())
        
        outfile = '{}.counts'.format(output)
        ref_range = None
    
    file_opened = False

    for ref_name in sorted(ref_names):
        header = ['ref_name', 'pos']
        if ref_to_seq:
            header.append('base')
        for sample in samples:
            header.append(sample)

        # for rrna case write to separate files for each ref_to_seq
        if ref_to_seq:
            outfile = '{}_{}.counts'.format(output, ref_name)
            ref_range = len(ref_to_seq[ref_name])

        if os.path.exists(outfile) and not file_opened:
            raise SystemExit('file {} already exists'.format(outfile))

        with open(outfile, 'a') as FHOUT:
            FHOUT.write(','.join(header) + '\n')    
            for pos in generate_rows(seen_pos[ref_name], ref_range):
                row = [ref_name, pos]
                # add the corresponding base is reference seq was given
                if ref_to_seq:
                    row.append(ref_to_seq[ref_name][pos])
                for sample in samples:
                    row.append(counts[sample][ref_name][pos])
                FHOUT.write(','.join([str(i) for i in row]) + '\n')
            if not ref_to_seq:
                file_opened = True

def main():
    parser = argparse.ArgumentParser(
        description= 'Count ends from bam alignments, requires samtools')
    parser.add_argument('bamfiles_samples', type=str, 
        help='file where each line: bam filename, sample name')
    parser.add_argument('output', type=str,
        help='output name (for example nm2_rRNA, or nm3_chrM)')
    parser.add_argument('-r', '--ref', type=str, 
        help='reference sequences in fasta format with headers that match bam references', 
        metavar='reference.fasta')
    parser.add_argument('-mm', '--mismatch',
        help='allow mismatches in bam alignments (default is no mismatches)',
        action='store_true')
    parser.add_argument('-f', '--samflag', type=int, 
        help='sam -f SAMFLAG (default is 131 for R2 reads)')
    parser.add_argument('--header',
        help='keep header in sam file (default is no header)',
        action='store_true')
    parser.add_argument('--keepsam',
        help='keep sam files (default is they are deleted)',
        action='store_true')
    parser.add_argument('--refnames', type=str,
        help='file containing <reference name>, <shortened name> (gi|12044..., 28S)')
    args = parser.parse_args()

    # default samflag is 131 = read paired, read mapped in proper pair, second in pair
    if not args.samflag:
        args.samflag = 131

    if args.ref:
        args.ref = fasta_to_dict(args.ref)
        # make sure ref seqs start with a hyphen so that they will be 0-based to match sam positions
        for header, seq in args.reg.items():
            if seq[0].upper() in 'ACGT':
                args.ref[header] = '-' + seq

    bamfiles_samples = read_bamfiles_samples(args.bamfiles_samples)

    samples = [sample for bamfile, sample in bamfiles_samples]

    samfiles_samples = bam_to_sam(bamfiles_samples, args.samflag, 
        no_mismatch=args.mismatch, keep_header=args.header)

    # count ends for each samfile
    counts = {}
    for samfile, sample in samfiles_samples:
        counts[sample] = count_ends(samfile, sample, args.ref)

    print(samples)
    print(counts.keys())
    print(counts['WT_BSA_input'])

    # write end counts to file
    write_counts(counts, args.output, samples, args.ref)

    # remove sam files

if __name__ == '__main__':
    main()