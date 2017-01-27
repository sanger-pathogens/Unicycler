"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module contains functions relating to Pilon polishing of a Unicycler assembly.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import os
import subprocess
import statistics
from collections import defaultdict
from .misc import load_fasta, reverse_complement, int_to_str, float_to_str, get_percentile_sorted


class CannotPolish(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


def polish_with_pilon(graph, bowtie2_path, bowtie2_build_path, pilon_path, java_path, samtools_path,
                      min_polish_size, polish_dir, verbosity, short_1, short_2, threads):
    """
    Runs Pilon on the graph to hopefully fix up small mistakes.
    """
    segments_to_polish = [x for x in graph.segments.values() if x.get_length() >= min_polish_size]
    if not segments_to_polish:
        raise CannotPolish('no segments are long enough to polish')

    polish_input_filename = os.path.join(polish_dir, 'polish.fasta')
    polish_fasta = open(polish_input_filename, 'w')
    for segment in segments_to_polish:
        polish_fasta.write('>' + str(segment.number) + '\n')
        polish_fasta.write(segment.forward_sequence)
        polish_fasta.write('\n')
    polish_fasta.close()

    # Prepare the FASTA for Bowtie2 alignment.
    bowtie2_build_command = [bowtie2_build_path, polish_input_filename, polish_input_filename]
    try:
        subprocess.check_output(bowtie2_build_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('bowtie2-build encountered an error:\n' + e.output.decode())
    if not any(x.endswith('.bt2') for x in os.listdir(polish_dir)):
        raise CannotPolish('bowtie2-build failed to build an index')

    # Perform the alignment with Bowtie2.
    raw_sam_filename = os.path.join(polish_dir, 'alignments_raw.sam')
    bowtie2_command = [bowtie2_path, '--end-to-end', '--very-sensitive', '--threads', str(threads),
                       '--no-discordant', '--no-mixed', '--no-unal', '-I', '0', '-X', '2000',
                       '-x', polish_input_filename, '-1', short_1, '-2', short_2,
                       '-S', raw_sam_filename]
    if verbosity > 0:
        print('Aligning short reads to assembly using Bowtie2')
    if verbosity > 1:
        print('  ' + ' '.join(bowtie2_command))
    try:
        subprocess.check_output(bowtie2_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('Bowtie2 encountered an error:\n' + e.output.decode())

    # CHECK TO MAKE SURE THAT alignments_raw.sam EXISTS AND LOOKS CORRECT.
    # if err:
    #     if verbosity > 1:
    #         print('\nAn error occurred during alignment:\n' + err.decode())
    #     raise CannotPolish('something')

    # Loop through alignments_raw.sam once to collect the insert sizes.
    insert_sizes = []
    raw_sam = open(raw_sam_filename, 'rt')
    for sam_line in raw_sam:
        try:
            insert_size = float(sam_line.split('\t')[8])
            if insert_size > 0.0:
                insert_sizes.append(insert_size)
        except (ValueError, IndexError):
            pass
    raw_sam.close()
    if not insert_sizes:
        raise CannotPolish('no read pairs aligned')
    insert_mean = statistics.mean(insert_sizes)
    if verbosity > 1:
        print('')
    if verbosity > 0:
        print('Mean fragment size =', float_to_str(insert_mean, 1) + ' bp')
    insert_sizes = sorted(insert_sizes)
    insert_size_5th = get_percentile_sorted(insert_sizes, 5.0)
    if verbosity > 1:
        print('Fragment size 5th percentile =', float_to_str(insert_size_5th, 0) + ' bp')
    insert_size_95th = get_percentile_sorted(insert_sizes, 95.0)
    if verbosity > 1:
        print('Fragment size 95th percentile =', float_to_str(insert_size_95th, 0) + ' bp')

    # Produce a new sam file including only the pairs with an appropriate insert size.
    filtered_sam_filename = os.path.join(polish_dir, 'alignments_filtered.sam')
    if verbosity > 0:
        print('Filtering alignments to fragment size range:',
              float_to_str(insert_size_5th, 0) + ' to ' + float_to_str(insert_size_95th, 0))
    filtered_sam = open(filtered_sam_filename, 'w')
    raw_sam = open(raw_sam_filename, 'rt')
    for sam_line in raw_sam:
        try:
            insert_size = abs(float(sam_line.split('\t')[8]))
            if insert_size_5th <= insert_size <= insert_size_95th:
                filtered_sam.write(sam_line)
        except (ValueError, IndexError):
            filtered_sam.write(sam_line)
    raw_sam.close()
    filtered_sam.close()

    # Sort the alignments.
    bam_filename = os.path.join(polish_dir, 'alignments.bam')
    samtools_sort_command = [samtools_path, 'sort', '-@', str(threads), '-o', bam_filename, '-O',
                             'bam', '-T', 'temp', filtered_sam_filename]
    if verbosity > 1:
        print('')
    if verbosity > 0:
        print('Sorting and indexing alignments')
    if verbosity > 1:
        print('  ' + ' '.join(samtools_sort_command))
    try:
        subprocess.check_output(samtools_sort_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('Samtools encountered an error:\n' + e.output.decode())

    # Index the alignments.
    samtools_index_command = [samtools_path, 'index', bam_filename]
    if verbosity > 1:
        print('  ' + ' '.join(samtools_index_command))
    try:
        subprocess.check_output(samtools_index_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('Samtools encountered an error:\n' + e.output.decode())

    # Polish with Pilon.
    if pilon_path.endswith('.jar'):
        pilon_command = [java_path, '-jar', pilon_path]
    else:
        pilon_command = [pilon_path]
    pilon_command += ['--genome', polish_input_filename, '--frags', bam_filename,
                      '--fix', 'bases', '--changes', '--outdir', polish_dir]
    if verbosity > 1:
        print('')
    if verbosity > 0:
        print('Running Pilon')
    if verbosity > 1:
        print('  ' + ' '.join(pilon_command))
    try:
        subprocess.check_output(pilon_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('Pilon encountered an error:\n' + e.output.decode())
    pilon_fasta_filename = os.path.join(polish_dir, 'pilon.fasta')
    pilon_changes_filename = os.path.join(polish_dir, 'pilon.changes')
    if not os.path.isfile(pilon_fasta_filename):
        raise CannotPolish('Pilon did not produce pilon.fasta')
    if not os.path.isfile(pilon_changes_filename):
        raise CannotPolish('Pilon did not produce pilon.changes')

    # Display Pilon changes.
    change_count = defaultdict(int)
    change_lines = defaultdict(list)
    total_count = 0
    pilon_changes = open(pilon_changes_filename, 'rt')
    for line in pilon_changes:
        try:
            seg_num = int(line.split(':')[0])
            change_count[seg_num] += 1
            total_count += 1
            change_lines[seg_num].append(line.strip())
        except ValueError:
            pass
    if verbosity > 0 and total_count == 0:
        print('No Pilon changes')
    elif verbosity == 1:
        print('Number of Pilon changes:', int_to_str(total_count))
    elif verbosity > 1:
        print('')
        seg_nums = sorted(graph.segments)
        polish_input_seg_nums = set(x.number for x in segments_to_polish)
        for seg_num in seg_nums:
            if seg_num in polish_input_seg_nums:
                count = change_count[seg_num]
                if count < 1:
                    continue
                print('Segment ' + str(seg_num) + ' (' +
                      int_to_str(graph.segments[seg_num].get_length()) + ' bp): ' +
                      int_to_str(count) + ' change' +
                      ('s' if count > 1 else ''))
                if verbosity > 2:
                    try:
                        changes = change_lines[seg_num]
                        changes = sorted(changes, key=lambda x:
                                         int(x.replace(' ', ':').replace('-', ':').split(':')[1]))
                        for change in changes:
                            print('  ' + change)
                    except (ValueError, IndexError):
                        pass

    # Replace segment sequences with Pilon-polished versions.
    pilon_results = load_fasta(pilon_fasta_filename)
    for header, sequence in pilon_results:
        if header.endswith('_pilon'):
            header = header[:-6]
        try:
            seg_num = int(header)
            segment = graph.segments[seg_num]
            segment.forward_sequence = sequence
            segment.reverse_sequence = reverse_complement(sequence)
        except (ValueError, KeyError):
            pass
