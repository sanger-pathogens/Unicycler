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
import shutil
from collections import defaultdict
from .misc import load_fasta, reverse_complement, int_to_str, underline, get_percentile_sorted, dim

from . import log


class CannotPolish(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


def polish_with_pilon_multiple_rounds(graph, bowtie2_path, bowtie2_build_path, pilon_path,
                                      java_path, samtools_path, min_polish_size, polish_dir,
                                      short_1, short_2, unpaired, threads, max_polish_count, keep):
    if not os.path.exists(polish_dir):
        os.makedirs(polish_dir)

    # To avoid issues with paths that contain spaces, we will move into the temporary directory
    # to run these commands.
    starting_dir = os.getcwd()
    os.chdir(polish_dir)

    insert_size_1st, insert_size_99th = \
        get_insert_size_range(graph, bowtie2_path, bowtie2_build_path, polish_dir, min_polish_size,
                              short_1, short_2, threads, keep)
    for i in range(max_polish_count):
        change_count = polish_with_pilon(graph, bowtie2_path, bowtie2_build_path, pilon_path,
                                         java_path, samtools_path, min_polish_size, polish_dir,
                                         short_1, short_2, unpaired, threads, insert_size_1st,
                                         insert_size_99th, i+1, keep)
        if change_count == 0:
            break

    if keep < 3 and os.path.exists(polish_dir):
        shutil.rmtree(polish_dir)
    os.chdir(starting_dir)


def get_insert_size_range(graph, bowtie2_path, bowtie2_build_path, polish_dir, min_polish_size,
                          short_1, short_2, threads, keep):
    """
    This function just does a quick alignment of the paired-end reads to figure out the 1st and
    99th percentiles for the insert size.
    """
    using_paired_reads = bool(short_1) and bool(short_2)
    if not using_paired_reads:
        return 0, 1000

    log.log('Aligning reads to find appropriate insert size range...')

    segments_to_polish = [x for x in graph.segments.values() if x.get_length() >= min_polish_size]
    if not segments_to_polish:
        raise CannotPolish('segments are too short')

    fasta_filename = '0_insert_size_check.fasta'
    sam_filename = '0_alignments.sam'

    with open(fasta_filename, 'w') as polish_fasta:
        for segment in segments_to_polish:
            polish_fasta.write('>' + str(segment.number) + '\n')
            polish_fasta.write(segment.forward_sequence)
            polish_fasta.write('\n')

    # Prepare the FASTA for Bowtie2 alignment.
    bowtie2_build_command = [bowtie2_build_path, fasta_filename, fasta_filename]
    log.log(dim('  ' + ' '.join(bowtie2_build_command)), 2)
    try:
        subprocess.check_output(bowtie2_build_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('bowtie2-build encountered an error:\n' + e.output.decode())
    if not any(x.endswith('.bt2') for x in os.listdir(polish_dir)):
        raise CannotPolish('bowtie2-build failed to build an index')

    # Perform the alignment with Bowtie2.
    bowtie2_command = [bowtie2_path, '-1', short_1, '-2', short_2, '-x', fasta_filename,
                       '--local', '--fast-local', '--threads', str(threads), '-I', '0',
                       '-X', '5000', '-S', sam_filename]
    log.log(dim('  ' + ' '.join(bowtie2_command)), 2)
    try:
        subprocess.check_output(bowtie2_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('Bowtie2 encountered an error:\n' + e.output.decode())

    insert_sizes = []
    with open(sam_filename, 'rt') as raw_sam:
        for sam_line in raw_sam:
            try:
                sam_parts = sam_line.split('\t')
                sam_flags = int(sam_parts[1])
                if sam_flags & 2:  # if read mapped in proper pair
                    insert_size = int(sam_parts[8])
                    if insert_size > 0.0:
                        insert_sizes.append(insert_size)
            except (ValueError, IndexError):
                pass
    if not insert_sizes:
        raise CannotPolish('no read pairs aligned')
    insert_sizes = sorted(insert_sizes)
    insert_size_1st = get_percentile_sorted(insert_sizes, 1.0)
    insert_size_99th = get_percentile_sorted(insert_sizes, 99.0)

    log.log('Insert size 1st percentile:  ' + str(insert_size_1st))
    log.log('Insert size 99th percentile: ' + str(insert_size_99th))
    log.log('')

    if keep < 3:
        for f in [fasta_filename, sam_filename, fasta_filename + '.1.bt2',
                  fasta_filename + '.2.bt2', fasta_filename + '.3.bt2', fasta_filename + '.4.bt2',
                  fasta_filename + '.rev.1.bt2', fasta_filename + '.rev.2.bt2']:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass

    return insert_size_1st, insert_size_99th


def polish_with_pilon(graph, bowtie2_path, bowtie2_build_path, pilon_path, java_path, samtools_path,
                      min_polish_size, polish_dir, short_1, short_2, unpaired, threads,
                      insert_size_1st, insert_size_99th, round_num, keep):
    """
    Runs Pilon on the graph to hopefully fix up small mistakes.
    """
    log.log(underline('Pilon polish round ' + str(round_num)))

    using_paired_reads = bool(short_1) and bool(short_2)
    using_unpaired_reads = bool(unpaired)

    input_filename = str(round_num) + '_polish_input.fasta'
    sam_filename = str(round_num) + '_alignments.sam'
    bam_filename = str(round_num) + '_alignments.bam'
    output_prefix = str(round_num) + '_pilon'
    pilon_fasta_filename = str(round_num) + '_pilon.fasta'
    pilon_changes_filename = str(round_num) + '_pilon.changes'
    pilon_output_filename = str(round_num) + '_pilon.out'

    segments_to_polish = [x for x in graph.segments.values() if x.get_length() >= min_polish_size]
    if not segments_to_polish:
        raise CannotPolish('no segments are long enough to polish')

    with open(input_filename, 'w') as polish_fasta:
        for segment in segments_to_polish:
            polish_fasta.write('>' + str(segment.number) + '\n')
            polish_fasta.write(segment.forward_sequence)
            polish_fasta.write('\n')

    # Prepare the FASTA for Bowtie2 alignment.
    bowtie2_build_command = [bowtie2_build_path, input_filename, input_filename]
    log.log(dim('  ' + ' '.join(bowtie2_build_command)), 2)
    try:
        subprocess.check_output(bowtie2_build_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('bowtie2-build encountered an error:\n' + e.output.decode())
    if not any(x.endswith('.bt2') for x in os.listdir(polish_dir)):
        raise CannotPolish('bowtie2-build failed to build an index')

    # Perform the alignment with Bowtie2.
    bowtie2_command = [bowtie2_path, '--local', '--very-sensitive-local', '--threads', str(threads),
                       '-I', str(insert_size_1st), '-X', str(insert_size_99th),
                       '-x', input_filename, '-S', sam_filename]
    if using_paired_reads:
        bowtie2_command += ['-1', short_1, '-2', short_2]
    if using_unpaired_reads:
        bowtie2_command += ['-U', unpaired]

    log.log(dim('  ' + ' '.join(bowtie2_command)), 2)
    try:
        subprocess.check_output(bowtie2_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('Bowtie2 encountered an error:\n' + e.output.decode())

    # Sort the alignments.
    samtools_sort_command = [samtools_path, 'sort', '-@', str(threads), '-o', bam_filename, '-O',
                             'bam', '-T', 'temp', sam_filename]
    log.log(dim('  ' + ' '.join(samtools_sort_command)), 2)
    try:
        subprocess.check_output(samtools_sort_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('Samtools encountered an error:\n' + e.output.decode())

    # Index the alignments.
    samtools_index_command = [samtools_path, 'index', bam_filename]
    log.log(dim('  ' + ' '.join(samtools_index_command)), 2)
    try:
        subprocess.check_output(samtools_index_command, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('Samtools encountered an error:\n' + e.output.decode())

    # Polish with Pilon.
    if pilon_path.endswith('.jar'):
        pilon_command = [java_path, '-jar', pilon_path]
    else:
        pilon_command = [pilon_path]
    pilon_command += ['--genome', input_filename, '--frags', bam_filename, '--changes',
                      '--output', output_prefix, '--outdir', polish_dir, '--fix', 'bases']
    log.log(dim('  ' + ' '.join(pilon_command)), 2)
    try:
        pilon_stdout = subprocess.check_output(pilon_command, stderr=subprocess.STDOUT)
        with open(pilon_output_filename, 'wb') as pilon_out:
            pilon_out.write(pilon_stdout)
    except subprocess.CalledProcessError as e:
        raise CannotPolish('Pilon encountered an error:\n' + e.output.decode())
    if not os.path.isfile(pilon_fasta_filename):
        raise CannotPolish('Pilon did not produce FASTA file')
    if not os.path.isfile(pilon_changes_filename):
        raise CannotPolish('Pilon did not produce changes file')

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
    log.log('', 2)
    if total_count == 0:
        log.log('No Pilon changes')
    else:
        log.log('Number of Pilon changes: ' + int_to_str(total_count))
        seg_nums = sorted(graph.segments)
        polish_input_seg_nums = set(x.number for x in segments_to_polish)
        for seg_num in seg_nums:
            if seg_num in polish_input_seg_nums:
                count = change_count[seg_num]
                if count < 1:
                    continue
                log.log('', 2)
                log.log('Segment ' + str(seg_num) + ' (' +
                        int_to_str(graph.segments[seg_num].get_length()) + ' bp): ' +
                        int_to_str(count) + ' change' +
                        ('s' if count > 1 else ''), 2)
                try:
                    changes = change_lines[seg_num]
                    changes = sorted(changes, key=lambda x:
                                     int(x.replace(' ', ':').replace('-', ':').split(':')[1]))
                    for change in changes:
                        log.log('  ' + dim(change), 2)
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

    log.log('', 2)

    if keep < 3:
        for f in [input_filename, sam_filename, bam_filename, bam_filename + '.bai',
                  pilon_fasta_filename, pilon_changes_filename, input_filename + '.1.bt2',
                  input_filename + '.2.bt2', input_filename + '.3.bt2',
                  input_filename + '.4.bt2', input_filename + '.rev.1.bt2',
                  input_filename + '.rev.2.bt2', pilon_output_filename]:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass

    return total_count
