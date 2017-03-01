#!/usr/bin/env python3
"""
This script runs minimap on Illumina contigs and long reads. It then produces two subsets of reads:
  1) Reads which overlap at least two single copy contigs. These are used for a miniasm assembly.
     These reads are also split in the output such that they never span a whole single copy contig.
  2) Reads which overlap at least one single copy contig. These are less useful for assembly but
     are still useful for Racon polishing.

Both subsets also include the Illumina contigs in FASTQ form.

Neither set includes:
  1) Reads which are contained within single copy contigs. These are not informative to assembly
     and possibly detrimental to Racon polishing.
  2) Reads which align to no single copy contigs. These are either contamination or reside in
     repeat regions. Either way - not useful for assembly or polishing.
"""

import argparse
import subprocess
import gzip


def main():
    args = get_arguments()

    reads = load_fastq(args.long_reads)
    print(str(len(reads)) + ' total reads')

    # Align the reads to the contigs with minimap.
    p = subprocess.Popen('minimap -w5 -L100 -m0 -t8 ' + args.contigs + ' ' + args.long_reads,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = p.communicate()
    paf_lines = stdout.decode().splitlines()

    # Load in the alignments and group them into read objects.
    aligned_reads = {}
    for paf_line in paf_lines:
        alignment = PafAlignment(paf_line)
        read_name = alignment.read_name
        if alignment.read_name not in aligned_reads:
            aligned_reads[alignment.read_name] = LongRead(alignment.read_name, alignment.read_length)
        aligned_reads[alignment.read_name].add_alignment(alignment)
    print(str(len(aligned_reads)) + ' reads aligned')

    one_or_more_count = 0
    with open('overlap_one_or_more.fastq', 'wt') as one_or_more:
        overlapping_one_or_more, overlapping_two_or_more = [], []
        for read in reads:
            read_name = read[0]
            if read_name in aligned_reads:
                aligned_read = aligned_reads[read_name]
                overlaps = aligned_read.overlapping_alignment_count()
                if overlaps >= 1:
                    one_or_more_count += 1
                    one_or_more.write('@')
                    one_or_more.write(read_name)
                    one_or_more.write('\n')
                    one_or_more.write(read[1])
                    one_or_more.write('\n+\n')
                    one_or_more.write(read[2])
                    one_or_more.write('\n')
    print(str(one_or_more_count) + ' overlap one or more contigs')

    two_or_more_count = 0
    with open('overlap_two_or_more.fastq', 'wt') as two_or_more:
        overlapping_one_or_more, overlapping_two_or_more = [], []
        for read in reads:
            read_name = read[0]
            if read_name in aligned_reads:
                aligned_read = aligned_reads[read_name]
                overlaps = aligned_read.overlapping_alignment_count()
                if overlaps >= 2:
                    two_or_more_count += 1
                    seq = read[1]
                    quals = read[2]
                    ranges = aligned_read.get_assembly_output_ranges()
                    for i, out_range in enumerate(ranges):
                        s = out_range[0]
                        e = out_range[1]
                        two_or_more.write('@')
                        two_or_more.write(read_name + '_' + str(i))
                        two_or_more.write('\n')
                        two_or_more.write(seq[s:e])
                        two_or_more.write('\n+\n')
                        two_or_more.write(quals[s:e])
                        two_or_more.write('\n')
    print(str(two_or_more_count) + ' overlap two or more contigs')









def get_arguments():
    parser = argparse.ArgumentParser(description='Subset long reads for miniasm bridging',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--contigs', type=str, required=True,
                        help='Reads shorter than this length (in bp) will not be outputted')
    parser.add_argument('--long_reads', type=str, required=True,
                        help='Reads with a mean qscore less than this will not be outputted')
    args = parser.parse_args()
    return args


class PafAlignment(object):
    """
    1   str  Query sequence name
    2   int  Query sequence length
    3   int  Query start coordinate (0-based)
    4   int  Query end coordinate (0-based)
    5   char `+' if query and target on the same strand; `-' if opposite
    6   str  Target sequence name
    7   int  Target sequence length
    8   int  Target start coordinate on the original strand
    9   int  Target end coordinate on the original strand
    10  int  Number of matching bases in the mapping
    11  int  Number bases, including gaps, in the mapping
    12  int  Mapping quality (0-255 with 255 for missing)
    """
    def __init__(self, paf_line):
        self.paf_line = paf_line
        line_parts = paf_line.split('\t')
        self.read_name = line_parts[0]
        self.read_length = int(line_parts[1])
        self.read_start = int(line_parts[2])
        self.read_end = int(line_parts[3])
        self.strand = line_parts[4]
        self.contig_name = line_parts[5]
        self.contig_length = int(line_parts[6])
        self.contig_start = int(line_parts[7])
        self.contig_end = int(line_parts[8])
        self.matching_bases = int(line_parts[9])
        self.bases = int(line_parts[10])

        self.read_end_gap = self.read_length - self.read_end
        self.contig_end_gap = self.contig_length - self.contig_end

    def __repr__(self):
        return str(self.read_start) + '-' + str(self.read_end) + ':' + \
            self.contig_name + ':' + str(self.contig_start) + '-' + str(self.contig_end) + \
            '(' + str(self.matching_bases) + ')'


class LongRead(object):
    def __init__(self, name, length):
        self.name = name
        self.length = length
        self.alignments = []

    def __repr__(self):
        return self.name + ' (' + str(self.length) + ' bp) ' + str(self.alignments) 

    def add_alignment(self, paf_alignment):
        self.alignments.append(paf_alignment)
        self.alignments = sorted(self.alignments, reverse=True, key=lambda x: x.matching_bases)
        best_matching_bases = self.alignments[0].matching_bases
        kept_alignments = []
        for alignment in self.alignments:
            if not alignments_overlap(alignment, kept_alignments) and alignment.matching_bases * 10 > best_matching_bases:
                kept_alignments.append(alignment)
        self.alignments = kept_alignments
        self.alignments = sorted(self.alignments, key=lambda x: x.read_start)

    def overlapping_alignment_count(self):
        """
        Counts the number of alignments for the read which overlap the end of a contig.
        """
        overlapping_alignments = 0
        for a in self.alignments:
            adjusted_contig_start = a.contig_start - a.read_start
            adjusted_contig_end = a.contig_end + a.read_end_gap
            if adjusted_contig_start < 0 or adjusted_contig_end >= a.contig_length:
                overlapping_alignments += 1
        return overlapping_alignments

    def get_assembly_output_ranges(self):
        """
        This function outputs the part(s) of the read which should be output as reads for assembly.
        """
        range_starts, range_ends = [], []
        for a in self.alignments:
            adjusted_contig_start = a.contig_start - a.read_start
            adjusted_contig_end = a.contig_end + a.read_end_gap
            contained_contig = adjusted_contig_start < 0 and adjusted_contig_end >= a.contig_length
            if contained_contig:
                range_starts.append(a.read_start + 100)
                range_ends.append(a.read_end - 100)
        range_starts = [0] + range_starts
        range_ends.append(self.length)
        return list(zip(range_starts, range_ends))


        if len(self.alignments) <= 1:
            return []
        for i in range(len(self.alignments) - 1):
            a1 = self.alignments[i]
            a2 = self.alignments[i+1]



def range_overlap(x1, x2, y1, y2):
    """
    Returns true if the range (x1, x2) overlaps with the range (y1, y2).
    """
    return x1 < y2 and y1 < x2


def alignments_overlap(a, other):
    return any(range_overlap(a.read_start, a.read_end, x.read_start, x.read_end) for x in other)


def load_fastq(filename):
    if get_compression_type(filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    reads = []
    with open_func(filename, 'rt') as fastq:
        for line in fastq:
            name = line.strip()[1:].split()[0]
            sequence = next(fastq).strip()
            _ = next(fastq)
            qualities = next(fastq).strip()
            reads.append((name, sequence, qualities))
    return reads


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        quit_with_error('cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        quit_with_error('cannot use zip format - use gzip instead')
    return compression_type


if __name__ == '__main__':
    main()
