#!/usr/bin/env python3
"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module contains the main script for the Unicycler assembly polisher. It is executed when a
user runs `unicycler_polish` (after installation) or `unicycler_polish-runner.py`.

It uses Pilon, Arrow, Racon and Freebayes (as appropriate for the input reads) to repeatedly polish
the assembly. The assembly quality is quantified using ALE and polishing will continue until no
more improvements are possible.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""


import argparse
import os
import math
import sys
from collections import defaultdict
from .misc import MyHelpFormatter, bold, quit_with_error, get_default_thread_count, \
    check_file_exists, int_to_str, float_to_str
from .minimap_alignment import load_minimap_alignments_basic
from .read_ref import load_long_reads
from . import log

try:
    from .cpp_wrappers import minimap_align_reads
except AttributeError as att_err:
    sys.exit('Error when importing C++ library: ' + str(att_err) + '\n'
             'Have you successfully built the library file using make?')


def main():
    full_command = ' '.join(('"' + x + '"' if ' ' in x else x) for x in sys.argv)
    args = get_arguments()
    print_intro_message(args, full_command)

    seq_dict, seq_names, _ = load_long_reads(args.input, silent=False,
                                             section_header='Loading sequences')
    log.log('')

    if args.trim > 0:
        log.log_section_header('Conducting alignments for trimming', single_newline=True)
        alignments = get_minimap_alignments_by_seq(args.input, args.reads, args.threads)

        log.log_section_header('Trimming sequences', single_newline=True)
        trim_sequences(seq_dict, seq_names, alignments, args.trim)




def get_arguments():
    parser = argparse.ArgumentParser(description='Unicycler-scrub - read trimming, chimera '
                                                 'detection and misassembly detection',
                                     formatter_class=MyHelpFormatter)

    parser.add_argument('-i', '--input', type=str, required=True,
                        help='These are the reads or assembly to be scrubbed (can be FASTA or '
                             'FASTQ format')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='The scrubbed reads or assembly will be saved to this file (will '
                             'have the same format as the --input file format)')
    parser.add_argument('-r', '--reads', type=str, required=True,
                        help='These are the reads used to scrub --input (can be FASTA or FASTQ '
                             'format, can be the same file as --input)')
    parser.add_argument('--trim', type=int, default=50,
                        help='The aggressiveness with which the input will be trimmed (0 to 100, '
                             'where 0 is no trimming and 100 is very aggressive trimming)')
    parser.add_argument('--split', type=int, default=50,
                        help='The aggressiveness with which the input will be split (0 to 100, '
                             'where 0 is no splitting and 100 is very aggressive splitting)')
    parser.add_argument('-t', '--threads', type=int, required=False,
                        default=get_default_thread_count(), help='Number of threads used')
    args = parser.parse_args()

    if args.trim < 0 or args.trim > 100:
        quit_with_error('--trim must be between 0 and 100 (inclusive)')
    if args.split < 0 or args.split > 100:
        quit_with_error('--split must be between 0 and 100 (inclusive)')
    if args.threads <= 0:
        quit_with_error('--threads must be at least 1')

    args.input = os.path.abspath(args.input)
    args.reads = os.path.abspath(args.reads)
    args.out = os.path.abspath(args.out)

    check_file_exists(args.input)
    check_file_exists(args.reads)

    return args


def print_intro_message(args, full_command):
    log.log_section_header('Starting Unicycler-scrub', single_newline=True)

    log.log_explanation('Unicycler-scrub uses local alignments (from minimap) to trim and/or '
                        'split an input sequence. Long reads can be used to scrub themselves (by '
                        'supplying the same reads to --input and --out) to trim the reads and split '
                        'chimeric reads. Long reads can also be used to check an assembly by '
                        'searching for misassembly points (by supplying the assembly to --input '
                        'and the long reads to --out).',
                        extra_empty_lines_after=0)
    log.log_explanation('For more information, please see https://github.com/rrwick/Unicycler')

    log.log('Command: ' + bold(full_command))
    log.log('')

    trim_level_str = '%3d' % args.trim
    if args.trim == 0.0:
        trim_level_str += ' (no trimming)'
    elif args.trim < 20.0:
        trim_level_str += ' (very conservative trimming)'
    elif args.trim < 40.0:
        trim_level_str += ' (conservative trimming)'
    elif args.trim < 60.0:
        trim_level_str += ' (moderate trimming)'
    elif args.trim < 80.0:
        trim_level_str += ' (aggressive trimming)'
    else:
        trim_level_str += ' (very aggressive trimming)'

    split_level_str = '%3d' % args.split
    if args.split == 0.0:
        split_level_str += ' (no splitting)'
    elif args.split < 20.0:
        split_level_str += ' (very conservative splitting)'
    elif args.split < 40.0:
        split_level_str += ' (conservative splitting)'
    elif args.split < 60.0:
        split_level_str += ' (moderate splitting)'
    elif args.split < 80.0:
        split_level_str += ' (aggressive splitting)'
    else:
        split_level_str += ' (very aggressive splitting)'

    log.log('Trim level:  ' + trim_level_str)
    log.log('Split level: ' +  split_level_str)


def get_minimap_alignments_by_seq(input, reads, threads):
    if input == reads:
        minimap_preset = 'scrub reads with reads'
    else:
        minimap_preset = 'scrub assembly with reads'
    minimap_alignments_str = minimap_align_reads(input, reads, threads, 3, minimap_preset)
    minimap_alignments = load_minimap_alignments_basic(minimap_alignments_str)
    log.log(str(len(minimap_alignments)) + ' alignments found')
    alignments_by_seq = defaultdict(list)
    for a in minimap_alignments:
        alignments_by_seq[a.ref_name].append(a)
    for seq_name in alignments_by_seq:
        alignments_by_seq[seq_name] = sorted(alignments_by_seq[seq_name], key=lambda x: x.ref_start)
    return alignments_by_seq


def trim_sequences(seq_dict, seq_names, alignments, trim_setting):
    length_before = 0
    for name in seq_names:
        length_before += seq_dict[name].get_length()
    log.log('Total length before trimming: ' + int_to_str(length_before) + ' bp')

    trim_setting = (trim_setting / 100.0) ** 2
    bases_trimmed = 0
    length_after = 0
    for name in seq_names:
        seq = seq_dict[name]
        seq_alignments = alignments[name]
        seq_length = seq.get_length()

        # Total up the read's mean depth, but we don't want to count alignments that are too local
        # (e.g. only part of the sequence is in common). So we scale each alignment's contribution
        # down by its overhang.
        mean_depth = 0.0
        for a in seq_alignments:
            fraction_ref = a.fraction_ref_aligned()
            overhang = a.get_total_overhang()
            alignment_size = a.num_bases
            scaling_factor = max(alignment_size - overhang, 0) / alignment_size
            fraction_ref *= scaling_factor
            mean_depth += fraction_ref

        target_depth = int(math.floor(trim_setting * mean_depth))
        trim_start_pos = 0
        trim_end_pos = seq_length
        if target_depth > 0:
            depth_changes = defaultdict(int)
            for a in seq_alignments:
                depth_changes[a.ref_start] += 1
                depth_changes[a.ref_end] -= 1
            depth_changes = sorted(depth_changes.items())

            # Find the first position which exceeds the target depth.
            depth = 0
            for pos, change in depth_changes:
                depth += change
                assert depth >= 0
                if depth >= target_depth:
                    trim_start_pos = pos
                    break

            # Go through the depth changes backwards to find the last position which exceeds the
            # target depth.
            depth = 0
            for pos, change in depth_changes[::-1]:
                depth -= change
                assert depth >= 0
                if depth >= target_depth:
                    trim_end_pos = pos
                    break

        seq.sequence = seq.sequence[trim_start_pos:trim_end_pos]
        seq.qualities = seq.qualities[trim_start_pos:trim_end_pos]
        bases_trimmed += trim_start_pos
        bases_trimmed += seq_length - trim_end_pos
        length_after += len(seq.sequence)

    log.log('Total bases trimmed:          ' +
            int_to_str(bases_trimmed, max_num=length_before) + ' bp')
    mean_bases_per_seq = bases_trimmed / len(seq_names)
    log.log('Mean bases trimmed per seq:   ' +
            float_to_str(mean_bases_per_seq, 1, max_num=length_before))
    log.log('Total length after trimming:  ' +
            int_to_str(length_after, max_num=length_before) + ' bp')
