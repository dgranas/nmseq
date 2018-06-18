#!/usr/bin/python

'''
Generates a bam_samples.txt file in current directory

Requires: 
argparse

Usage:
python CreateBamSample.py -d [bam_directory]

Optional arguments:
-d [bam_directory]: directory containing bam files
default is to look in the current directory

Each bam file is given a sample name
Sample names can be given in 2 ways

1. Break up the bam filenames with a delimiter

    User enters the delimiter and the fields they want to keep
    NmSeq2_WT_BSA_input_filtered.bam with a delimiter of '_' and fields 2,3,4
    gives a sample name of WT_BSA_input
    This will be applied to all bam files

2. User names each bam file individually

The csv file bam_samples.txt is created where each line has
[bam filename], [sample name]

This can then be used as input for CountEnds.py
'''

import os
import argparse
from builtins import input # to make py2 and py3 input compatible

def find_bam_files(bam_dir = '.'):
    '''
    Looks for bam files in a directory
    Gets the sample names for each file
    This is done either by parsing with a delimiter
    Or by having the user enter the name for each file
    Returns list of tuples [(bamfile, sample)]
    '''
    # make sure there isn't already a bam_sample file in current directory
    if os.path.exists('bam_sample.txt'):
        raise SystemExit('There is already a bam_sample.txt file in the directory')
    bam_samples = []
    parse_asked = False # once we ask about parsing, set to True so don't ask again
    parse_flag = False # set to True if we are parsing all the files

    for filename in os.listdir(bam_dir):
        if filename.endswith('.bam'):
            print(filename)
            if not parse_asked:
                parse_asked = True

                delim = input('If you want to name the sample by breaking up the filename using a delimiter,\nenter the delimiter to use (or press Enter to skip this): ')
                if delim:
                    print('Enter the fields you want to include as comma-separated numbers, starting with 1')
                    print('For example, enter 2,3,4 to get WT_BSA_input from NmSeq_WT_BSA_input_filtered.bam')
                    
                    fields = []
                    while not fields:
                        fields = [int(i) for i in input('Enter fields to use: ').split(',')]
                    parse_flag = True

            if parse_flag:
                sample = '_'.join([filename.split(delim)[i-1] for i in fields])
            else:
                sample = input('Enter sample name (or press Enter to exclude): ')
            if sample:
                bam_samples.append((os.path.join(bam_dir, filename), sample))

    return bam_samples

def write_bam_samples(bam_samples):
    '''
    Get the bam files and their sample names from user input
    Writes the file bam_sample.txt
    '''
    bam_samples = sorted(bam_samples) # sort by filename

    print('\nEnter the desired order for the samples using comma-separated positions')
    print('For example entering 3,4,1,2 for samples A,B,C,D would get the order C,D,A,B')
    print('Samples: ')
    for i, bam_sample in enumerate(bam_samples):
        print('{}\t{}'.format(i+1, bam_sample[1]))
    order = input('Desired order (or press Enter to keep this order): ')
    if order:
        bam_samples = [bam_samples[int(i)-1] for i in order.split(',')]
        print('Samples reordered:')
        for i, bam_sample in enumerate(bam_samples):
            print('{}\t{}'.format(i+1, bam_sample[1]))

    # write out the new bam_sample.txt file
    with open('bam_sample.txt', 'w') as f:
        for bamfile, sample in bam_samples:
            f.write('{}, {}\n'.format(bamfile, sample))

    print('The file bam_sample.txt was created in current directory')

def main():
    parser = argparse.ArgumentParser(
        description= 'Generates a bam_samples.txt file')
    parser.add_argument('-d', '--directory', type=str,
        help='directory containing bam files')
    args = parser.parse_args()

    if not args.directory:
        args.directory = '.'

    bam_samples = find_bam_files(args.directory)

    # case were no bam files were found
    if not bam_samples:
        raise SystemExit('No bam files found in directory {}'.format(args.directory))

    write_bam_samples(bam_samples)

if __name__ == '__main__':
    main()
