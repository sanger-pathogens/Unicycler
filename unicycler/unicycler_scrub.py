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
import subprocess
from collections import defaultdict
from .misc import MyHelpFormatter, bold, quit_with_error, get_default_thread_count, \
    check_file_exists, int_to_str, float_to_str, get_sequence_file_type
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
    input_type = get_sequence_file_type(args.input)

    seq_dict, seq_names, _ = load_long_reads(args.input, silent=False,
                                             section_header='Loading sequences')
    log.log('')

    log.log_section_header('Conducting alignments', single_newline=True)
    alignments = get_minimap_alignments_by_seq(args.input, args.reads, args.threads)

    # Trim the sequences based on
    if args.trim > 0:
        log.log_section_header('Trimming sequences', single_newline=True)
        trim_sequences(seq_dict, seq_names, alignments, args.trim)
    else:
        for seq in seq_dict.values():
            seq.trim_start_pos = 0
            seq.trim_end_pos = seq.get_length()

    # Split the sequences where they look like chimeras/misassemblies.
    if args.split > 0:
        log.log_section_header('Splitting sequences', single_newline=True)
        split_sequences(seq_dict, seq_names, alignments, args.split, args.min_split_size,
                        args.discard_chimeras)
    else:
        for seq in seq_dict.values():
            seq.positive_score_ranges = [(0, seq.get_length())]

    # Combine the trimming with the splitting to get final output ranges for each sequence.
    for name in seq_names:
        seq = seq_dict[name]
        seq.final_ranges = []
        for s, e in seq.positive_score_ranges:
            s = max(s, seq.trim_start_pos)
            e = min(e, seq.trim_end_pos)
            if e - s >= args.min_split_size:
                seq.final_ranges.append((s, e))

    output_sequences(args.out, seq_names, seq_dict, input_type)


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
    parser.add_argument('--min_split_size', type=int, default=1000,
                        help='Parts of split sequences will only be outputted if they are at '
                             'least this big')
    parser.add_argument('--discard_chimeras', action='store_true',
                        help='If set, chimeric reads will be discarded instead of split')
    parser.add_argument('-t', '--threads', type=int, required=False,
                        default=get_default_thread_count(), help='Number of threads used')
    args = parser.parse_args()

    if args.trim < 0 or args.trim > 100:
        quit_with_error('--trim must be between 0 and 100 (inclusive)')
    if args.split < 0 or args.split > 100:
        quit_with_error('--split must be between 0 and 100 (inclusive)')
    if args.trim == 0 and args.split == 0:
        quit_with_error('--trim and --split cannot both be 0 (there would be nothing left to do)')
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
    log.log(int_to_str(len(minimap_alignments)) + ' alignments found')
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
        mean_depth = get_mean_seq_depth(seq_alignments)
        target_depth = int(math.floor(trim_setting * mean_depth))

        seq.trim_start_pos = 0
        seq.trim_end_pos = seq_length

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
                    seq.trim_start_pos = pos
                    break

            # Go through the depth changes backwards to find the last position which exceeds the
            # target depth.
            depth = 0
            for pos, change in depth_changes[::-1]:
                depth -= change
                assert depth >= 0
                if depth >= target_depth:
                    seq.trim_end_pos = pos
                    break

        bases_trimmed += seq.trim_start_pos
        bases_trimmed += seq_length - seq.trim_end_pos
        length_after += seq.trim_end_pos - seq.trim_start_pos

    log.log('Total bases trimmed:          ' +
            int_to_str(bases_trimmed, max_num=length_before) + ' bp')
    mean_bases_per_seq = bases_trimmed / len(seq_names)
    log.log('Mean bases trimmed per seq:   ' +
            float_to_str(mean_bases_per_seq, 1, max_num=length_before))
    log.log('Total length after trimming:  ' +
            int_to_str(length_after, max_num=length_before) + ' bp')


def split_sequences(seq_dict, seq_names, alignments, split_setting, min_split_size,
                    discard_chimeras):
    split_setting /= 100.0
    max_relevant_distance = 1000

    for name in seq_names:
        seq = seq_dict[name]
        seq_alignments = alignments[name]
        seq_length = seq.get_length()

        # Each position in the sequence is scored based on the alignments around it.
        starting_score = 1.0 - split_setting
        scores = [starting_score] * seq_length
        start_overhang_scores = [0.0] * seq_length
        end_overhang_scores = [0.0] * seq_length

        for a in seq_alignments:

            # Alignments which span the position earn it positive points (because they suggest
            # contiguity). The larger the span (up to a maximum), the more points are earned.
            a_start, a_end = a.ref_start, a.ref_end
            for pos in range(a_start, a_end):
                distance_from_end = min(pos - a_start, a_end - pos)
                relevant_distance = min(distance_from_end, max_relevant_distance)
                score = relevant_distance / max_relevant_distance
                score /= split_setting  # lower split settings give larger positive scores
                scores[pos] += score

            # Tally up the score for start overlaps. Scores are larger for positions in larger
            # overlaps and positions close to where the overlap begins.
            start_overhang_size = min(a.get_start_overhang(), max_relevant_distance)
            start_overhang_rel_size = start_overhang_size / max_relevant_distance
            for pos in range(a_start - start_overhang_size, a_start):
                distance_from_clip = a_start - pos
                score = (start_overhang_size - distance_from_clip) / start_overhang_size
                score *= start_overhang_rel_size
                start_overhang_scores[pos] -= score

            # Do the same for end overlaps.
            end_overhang_size = min(a.get_end_overhang(), max_relevant_distance)
            end_overhang_rel_size = end_overhang_size / max_relevant_distance
            for pos in range(a_end, a_end + end_overhang_size):
                distance_from_clip = pos - a_end
                score = (end_overhang_size - distance_from_clip) / end_overhang_size
                score *= end_overhang_rel_size
                end_overhang_scores[pos] -= score

        # The final score comes from the positive contribution of spanning alignments and the
        # negative contribution of overhangs. But there needs to be both a start overhang and
        # end overhang to really matter.
        seq.positive_score_ranges = []
        positive_range_start = 0
        in_positive_range = True
        for pos in range(seq_length):
            overhang_score = max(start_overhang_scores[pos], end_overhang_scores[pos])
            scores[pos] += overhang_score
            if scores[pos] >= 0.0:
                if not in_positive_range:
                    in_positive_range = True
                    positive_range_start = pos
            else:
                if in_positive_range:
                    in_positive_range = False
                    seq.positive_score_ranges.append((positive_range_start, pos))
        if in_positive_range:
            seq.positive_score_ranges.append((positive_range_start, seq_length))

        is_chimera = seq.positive_score_ranges != [(0, seq_length)]

        if is_chimera:
            if discard_chimeras:
                seq.positive_score_ranges = []
            else:
                # Clean up small ranges.
                seq.positive_score_ranges = [x for x in seq.positive_score_ranges
                                             if x[1] - x[0] >= min_split_size]

            if seq.positive_score_ranges:
                log.log('\nSplit ' + seq.name + ' into pieces: ' +
                        get_read_range_str(seq.positive_score_ranges))
                # print(seq.sequence)  # TEMP
            else:
                log.log('\nDiscarded ' + seq.name)
        else:
            log.log('.', end='')


def get_mean_seq_depth(seq_alignments):
    """
    Total up the read's mean depth, but we don't want to count alignments that are too local (e.g.
    only part of the sequence is in common). So we scale each alignment's contribution down by its
    overhang.
    """
    mean_depth = 0.0
    for a in seq_alignments:
        fraction_ref = a.fraction_ref_aligned()
        overhang = a.get_total_overhang()
        alignment_size = a.num_bases
        scaling_factor = max(alignment_size - overhang, 0) / alignment_size
        fraction_ref *= scaling_factor
        mean_depth += fraction_ref
    return mean_depth


def get_read_range_str(ranges):
    return ', '.join([str(x[0]) + '-' + str(x[1]) for x in ranges])


def output_sequences(output, seq_names, seq_dict, out_format):
    log.log('\n')
    gzipped_out = output.endswith('.gz')
    if gzipped_out:
        out_filename = 'TEMP_' + str(os.getpid()) + '.fastq'
    else:
        out_filename = output
    with open(out_filename, 'wt') as out:
        for name in seq_names:
            seq = seq_dict[name]
            include_piece_number = len(seq.final_ranges) > 1
            for i, range in enumerate(seq.final_ranges):
                if out_format == 'FASTA':
                    out_str = get_fasta(seq.name, range, seq.sequence, i, include_piece_number)
                else:  # FASTQ
                    out_str = get_fastq(seq.name, range, seq.sequence, seq.qualities, i,
                                        include_piece_number)
                out.write(out_str)
    if gzipped_out:
        subprocess.check_output('gzip -c ' + out_filename + ' > ' + output,
                                stderr=subprocess.STDOUT, shell=True)
        os.remove(out_filename)
    log.log('Saved scrubbed sequences to ' + os.path.abspath(output))


def get_fasta(name, range, sequence, i, include_piece_number):
    s, e = range
    parts = ['>', name]
    if include_piece_number:
        parts.append('_' + str(i+1))
    parts += ['\n', sequence[s:e], '\n']
    return ''.join(parts)


def get_fastq(name, range, sequence, qualities, i, include_piece_number):
    s, e = range
    parts = ['@', name]
    if include_piece_number:
        parts.append('_' + str(i+1))
    parts += ['\n', sequence[s:e], '\n+\n', qualities[s:e], '\n']
    return ''.join(parts)
