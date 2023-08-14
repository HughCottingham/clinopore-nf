#!/usr/bin/env python3
import argparse
import pathlib
import math
import statistics


from Bio.SeqIO.FastaIO import SimpleFastaParser


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--assembly_fps', required=True, type=pathlib.Path,
            help='Input FASTA file')
    parser.add_argument('-i', '--id', required=False, type=str, help='ID to put in assembly column, only used if only one sample')

    # Ensure that input file exists
    args = parser.parse_args()
    if not args.assembly_fps.exists():
        parser.error('Input file %s does not exist' % assembly_fp)

    return args


def main():
    # Get commandline arguments
    args = get_arguments()

    # Print header
    header = ['assembly', 'contig_number', 'n50', 'q1', 'q2', 'q3', 'mean', 'smallest', 'largest', 'length','circular_contigs','circular_chromosome']
    print(*header, sep='\t')
    # Print stats as we iterate
    ordered_stats = get_assembly_stats(args.assembly_fps, args.id)
    print(*ordered_stats, sep='\t')


def get_assembly_stats(assembly_fp, id):
    # Get contig lengths
    with assembly_fp.open('r') as f:
        contig_lengths = [len(s) for d, s in SimpleFastaParser(f)]

    # Calculate stats
    contig_number = len(contig_lengths)
    length = sum(contig_lengths)
    smallest = min(contig_lengths)
    largest = max(contig_lengths)
    mean = int(round(statistics.mean(contig_lengths), 0))
    q1, q2, q3 = calculate_quartiles(contig_lengths)
    # If we only have one contig, use that as n50 (otherwise n50 calc fails)
    if contig_number == 1:
        n50 = largest
    else:
        n50 = calculate_n50(contig_lengths, length/2)
    with assembly_fp.open('r') as f:
        circular_contigs=[]
        circular_chromosomes=[]
        for line in f:
            line=str(line)
            if '>' in line and 'circular=Y' in line and 'contig_1' in line:
                circular_contigs.append(line)
                circular_chromosomes.append(line)
            elif '>' in line and 'circular=Y' in line:
                circular_contigs.append(line)
        num_circular_contigs=len(circular_contigs)
        num_circular_chromosomes=len(circular_chromosomes)
        if largest>1000000 and num_circular_chromosomes>0:
            circular_chromosome='True'
        else:
            circular_chromosome='False'

    # Return ordered stats
    return id, contig_number, n50, q1, q2, q3, mean, smallest, largest, length, num_circular_contigs, circular_chromosome


def calculate_quartiles(lengths):
    # Set up
    np = [0.25, 0.50, 0.75]
    x = sorted(lengths)
    n = len(x)
    # Get bounds
    indices = [(n - 1) * p for p in np]
    lo = [math.floor(i) for i in indices]
    hi = [math.ceil(i) for i in indices]
    qs = [x[i] for i in lo]
    # Update if required and then return
    for i in range(len(indices)):
        if not indices[i] > lo[i]:
            continue
        h = indices[i] - lo[i]
        qs[i] = (1 - h) * qs[i] + h * x[hi[i]]
    return qs


def calculate_n50(lengths, median):
    csum = 0
    for length in sorted(lengths, reverse=True):
        # If cumulative sum exceeds median, return last length
        if csum > median:
            return prev_length
        # Cumulative sum
        csum += length
        # Set current to prev_length and then loop
        prev_length = length


if __name__ == '__main__':
    main()