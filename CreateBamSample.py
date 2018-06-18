#!/usr/bin/python

'''
Generates a bam_samples.txt file
Optional argument is to supply the folder with the bam files
Default is to look in current directory
'''

import os
import argparse

def find_bam_files(bam_dir = '.'):
    '''
    Looks for bam files in a directory
    Gets the sample names for each file
    This is done either by parsing with a delimiter
    Or by having the user enter the name for each file
    Returns list of tuples [(bamfile, sample)]
    '''
    # make sure there isn't already a bam_sample file in this directory
    if os.path.exists(os.path.join(bam_dir, 'bam_sample.txt')):
        raise SystemExit('There is already a bam_sample.txt file in the directory')
    bam_samples = []
    parse_asked = False # once we ask about parsing, set to True so don't ask again
    parse_flag = False # set to True if we are parsing all the files

    for filename in os.listdir(bam_dir):
        if filename.endswith('.bam'):
            print(filename)
            if not parse_asked:
                parse_asked = True

                delim = raw_input('If you want to name the sample by breaking up \
                    the filename using a delimiter,\n enter the delimiter to use \
                    (or press Enter to skip this): ')
                if delim:
                    print('Enter the fields you want to include as comma-separated \
                        numbers, starting with 1')
                    print('For example, enter 2,3,4 to get WT_BSA_input from \
                        NmSeq_WT_BSA_input_filtered.bam')
                    
                    fields = []
                    while not fields:
                        fields = [int(i) for i in raw_input('Enter fields to use: ').split(',')]
                    parse_flag = True

            if parse_flag:
                sample = '_'.join([filename.split(delim)[i-1] for i in fields])
            else:
                sample = raw_input('Enter sample name (or press Enter to exclude): ')
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
    order = raw_input('Desired order (or press Enter to keep this order): ')
    if order:
        bam_samples = [bam_samples[int(i)-1] for i in order.split(',')]
        print('Samples reordered:')
        for i, bam_sample in enumerate(samples):
            print('{}\t{}'.format(i+1, bam_sample[1]))

    # write out the new bam_sample.txt file
    with open('bam_sample.txt', 'w') as f:
        for bamfile, sample in bam_samples:
            f.write('{}, {}\n'.format(bamfile, sample))
            break

    print('The file bam_sample.txt was created')

def main():
    parser = argparse.ArgumentParser(
        description= 'Generates a bam_samples.txt file')
    parser.add_argument('-d', '--directory', type=str,
        help='directory containing bam files')
    args = parser.parse_args()

    if not args.directory:
        args.directory = '.'

    bam_samples = find_bam_files(args.directory)

    # if nothing found, ask user for the directory path
    if not bam_samples:
        raise SystemExit('No bam files found in directory {}'.format(bam_dir))

    write_bam_samples(bam_samples)

if __name__ == '__main__':
    main()
    