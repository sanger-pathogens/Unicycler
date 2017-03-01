"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module contains functions relating to SPAdes assembly.

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
import gzip
import shutil
import statistics
from .misc import round_to_nearest_odd, get_compression_type, int_to_str, quit_with_error,\
    strip_read_extensions, bold, dim, print_table, get_left_arrow
from .assembly_graph import AssemblyGraph
from . import log


def get_best_spades_graph(short1, short2, short_unpaired, out_dir, read_depth_filter, verbosity,
                          spades_path, threads, keep, kmer_count, min_k_frac, max_k_frac,
                          no_spades_correct, expected_linear_seqs):
    """
    This function tries a SPAdes assembly at different k-mers and returns the best.
    'The best' is defined as the smallest dead-end count after low-depth filtering.  If multiple
    graphs have the same dead-end count (e.g. zero!) then the highest kmer is used.
    """
    spades_dir = os.path.join(out_dir, 'spades_assembly')
    if not os.path.exists(spades_dir):
        os.makedirs(spades_dir)

    if no_spades_correct:
        reads = (short1, short2, short_unpaired)
    else:
        reads = spades_read_correction(short1, short2, short_unpaired, spades_dir, threads,
                                       spades_path)
    kmer_range = get_kmer_range(short1, short2, spades_dir, kmer_count, min_k_frac, max_k_frac)
    assem_dir = os.path.join(spades_dir, 'assembly')
    log.log_section_header('SPAdes assemblies')

    # Conduct a SPAdes assembly for each k-mer and score them to choose
    # the best.
    if verbosity > 1:
        spades_results_table = [['K-mer', 'Segments', 'Links', 'Total length', 'N50',
                                 'Longest segment', 'Dead ends', 'Score']]
    else:
        spades_results_table = [['K-mer', 'Segments', 'Dead ends', 'Score']]
    best_score = 0.0
    best_kmer = 0
    best_graph_filename = ''

    graph_files, insert_size_mean, insert_size_deviation = \
        spades_assembly(reads, assem_dir, kmer_range, threads, spades_path)

    existing_graph_files = [x for x in graph_files if x is not None]
    median_segment_count = statistics.median(count_segments_in_spades_fastg(x)
                                             for x in existing_graph_files)

    for graph_file, kmer in zip(graph_files, kmer_range):
        table_line = [int_to_str(kmer)]

        if graph_file is None:
            table_line += [''] * (7 if verbosity > 1 else 2)
            table_line.append('failed')
            spades_results_table.append(table_line)
            continue

        assembly_graph = AssemblyGraph(graph_file, kmer, paths_file=None,
                                       insert_size_mean=insert_size_mean,
                                       insert_size_deviation=insert_size_deviation)

        # If this graph has way too many segments, then we will just skip it because very complex
        # graphs take forever to clean up.
        # TO DO: I can remove this awkward hack if I make the graph cleaning more efficient.
        if len(assembly_graph.segments) > 4 * median_segment_count:
            table_line += [''] * (6 if verbosity > 1 else 2)
            table_line.append('too complex')
            spades_results_table.append(table_line)
            continue

        assembly_graph.clean(read_depth_filter)
        clean_graph_filename = os.path.join(spades_dir, ('k%03d' % kmer) + '_assembly_graph.gfa')
        assembly_graph.save_to_gfa(clean_graph_filename, verbosity=2)

        segment_count = len(assembly_graph.segments)
        dead_ends = assembly_graph.total_dead_end_count()

        # If the user is expecting some linear sequences, then the dead end count can be adjusted
        # down so expected dead ends don't penalise this k-mer.
        adjusted_dead_ends = max(0, dead_ends - (2 * expected_linear_seqs))
        if segment_count == 0:
            score = 0.0
        else:
            score = 1.0 / (segment_count * ((adjusted_dead_ends + 1) ** 2))

        # Prepare the table line for this k-mer graph.
        table_line += [int_to_str(segment_count)]
        if verbosity > 1:
            n50, shortest, lower_quartile, median, upper_quartile, longest = \
                assembly_graph.get_contig_stats()
            table_line += [int_to_str(assembly_graph.get_total_link_count()),
                           int_to_str(assembly_graph.get_total_length()),
                           int_to_str(n50), int_to_str(longest)]
        table_line += [int_to_str(dead_ends), '{:.2e}'.format(score)]
        spades_results_table.append(table_line)

        if score > best_score:
            best_kmer = kmer
            best_score = score
            best_graph_filename = graph_file

    log.log('', 2)

    if not best_kmer:
        quit_with_error('none of the SPAdes graphs were suitable for scaffolding in Unicycler')

    # Load and clean the best graph.
    assembly_graph = AssemblyGraph(best_graph_filename, best_kmer, paths_file=None,
                                   insert_size_mean=insert_size_mean,
                                   insert_size_deviation=insert_size_deviation)
    assembly_graph.clean(read_depth_filter)
    clean_graph_filename = os.path.join(spades_dir, 'k' + str(best_kmer) + '_assembly_graph.gfa')
    assembly_graph.save_to_gfa(clean_graph_filename, verbosity=2)

    if best_score == 0.0:
        quit_with_error('none of the SPAdes assemblies produced assembled sequence')

    # Print the SPAdes result table, highlighting the best k-mer in green.
    log.log_section_header('SPAdes assembly graph summary', 2)
    best_kmer_row = [x[0] for x in spades_results_table].index(int_to_str(best_kmer))
    print_table(spades_results_table, alignments='RRRRRRRR', indent=0,
                row_colour={best_kmer_row: 'green'},
                row_extra_text={best_kmer_row: ' ' + get_left_arrow() + 'best'})

    # Clean up.
    if keep < 3 and os.path.isdir(spades_dir):
        log.log('\nDeleting ' + spades_dir + '/')
        shutil.rmtree(spades_dir)

    return assembly_graph


def spades_read_correction(short1, short2, short_unpaired, spades_dir, threads, spades_path):
    """
    This runs SPAdes with the --only-error-correction option.
    """
    log.log_section_header('SPAdes read error correction')

    # If the corrected reads already exist, then we just use them and proceed.
    corrected_1 = os.path.join(spades_dir, 'corrected_1.fastq.gz')
    corrected_2 = os.path.join(spades_dir, 'corrected_2.fastq.gz')
    corrected_u = os.path.join(spades_dir, 'corrected_u.fastq.gz')
    if os.path.isfile(corrected_1) and os.path.isfile(corrected_2):
        log.log('Corrected reads already exist. Will use these reads instead of running SPAdes '
              'error correction:')
        log.log('  ' + corrected_1)
        log.log('  ' + corrected_2)
        if os.path.isfile(corrected_u):
            log.log('  ' + corrected_u)
        return corrected_1, corrected_2, corrected_u

    # If the corrected reads don't exist, then we run SPAdes in error correction only mode.
    read_correction_dir = os.path.join(spades_dir, 'read_correction')
    command = [spades_path, '-1', short1, '-2', short2]
    if short_unpaired:
        command += ['-s', short_unpaired]
    command += ['-o', read_correction_dir, '--threads', str(threads), '--only-error-correction']
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    while process.poll() is None:
        spades_output = process.stdout.readline().rstrip().decode()
        if 'Command line:' in spades_output:
            spades_output = ' '.join(spades_output.split()).replace('Command line: ', '')
            log.log('Command: ' + bold(spades_output))
            log.log('', 2)
        elif spades_output:
            log.log(dim(spades_output), 2)

    spades_error = process.stderr.readline().strip().decode()
    if spades_error:
        quit_with_error('SPAdes encountered an error: ' + spades_error)

    # Read error correction should be done now, so copy the correct read files to a more permanent
    # location.
    short1_no_extension = strip_read_extensions(short1)
    short2_no_extension = strip_read_extensions(short2)
    corrected_dir = os.path.join(read_correction_dir, 'corrected')
    files = os.listdir(corrected_dir)
    for spades_file in files:
        file_path = os.path.join(corrected_dir, spades_file)
        if short1_no_extension in spades_file:
            shutil.move(file_path, corrected_1)
        elif short2_no_extension in spades_file:
            shutil.move(file_path, corrected_2)
        elif '_unpaired' in spades_file:
            shutil.move(file_path, corrected_u)
    shutil.rmtree(read_correction_dir)

    if not os.path.isfile(corrected_1) or not os.path.isfile(corrected_2):
        quit_with_error('SPAdes read error correction failed')

    if not os.path.isfile(corrected_u):
        corrected_u = None

    log.log('')
    log.log('Corrected reads:')
    log.log('  ' + corrected_1)
    log.log('  ' + corrected_2)
    if corrected_u:
        log.log('  ' + corrected_u)

    return corrected_1, corrected_2, corrected_u


def spades_assembly(read_files, out_dir, kmers, threads, spades_path, just_last=False):
    """
    This runs a SPAdes assembly, possibly continuing from a previous assembly.
    """
    short1 = read_files[0]
    short2 = read_files[1]
    unpaired = read_files[2]
    kmer_string = ','.join([str(x) for x in kmers])
    command = [spades_path, '-o', out_dir, '-k', kmer_string, '--threads', str(threads)]
    if just_last:
        command += ['--restart-from', 'k' + str(kmers[-1])]
    else:
        command += ['--only-assembler', '-1', short1, '-2', short2]
        if unpaired:
            command += ['-s', unpaired]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    insert_size_mean = None
    insert_size_deviation = None
    while process.poll() is None:
        spades_output = process.stdout.readline().rstrip().decode()
        if spades_output:
            # Some SPAdes output lines use tabs where spaces would look better. Fix those up here
            # for aesthetics.
            if spades_output.startswith('Command line:') or \
                    spades_output.startswith('Restored from Command line:'):
                spades_output = ' '.join(spades_output.split())

            if spades_output.startswith('Command line:'):
                spades_output = spades_output.replace('Command line: ', '')
                log.log('Command: ' + bold(spades_output), 2)
                log.log('', 2)
            elif 'Running assembler: K' in spades_output:
                log.log(spades_output, 2)
            elif spades_output:
                log.log(dim(spades_output), 2)
        if 'Insert size =' in spades_output and 'deviation = ' in spades_output:
            try:
                insert_size_mean = float(spades_output.split('Insert size = ')[-1].split(',')[0])
                insert_size_deviation = float(spades_output.split('deviation = ')[-1].split(',')[0])
            except ValueError:
                pass
    log.log('', 2)

    spades_error = process.stderr.readline().strip().decode()
    if spades_error:
        quit_with_error('SPAdes encountered an error: ' + spades_error)

    if just_last:
        graph_file = os.path.join(out_dir, 'K' + str(kmers[-1]), 'assembly_graph.fastg')
        return graph_file, insert_size_mean, insert_size_deviation
    else:
        graph_files = []
        for kmer in kmers:
            graph_file = os.path.join(out_dir, 'K' + str(kmer), 'assembly_graph.fastg')
            if os.path.isfile(graph_file):
                parent_dir = os.path.dirname(out_dir)
                copied_graph_file = os.path.join(parent_dir,
                                                 ('k%03d' % kmer) + '_assembly_graph.fastg')
                shutil.copyfile(graph_file, copied_graph_file)
                graph_files.append(copied_graph_file)
            else:
                graph_files.append(None)
        return graph_files, insert_size_mean, insert_size_deviation


def get_kmer_range(reads_1_filename, reads_2_filename, spades_dir, kmer_count, min_kmer_frac,
                   max_kmer_frac):
    """
    Uses the read lengths to determine the k-mer range to be used in the SPAdes assembly.
    """
    log.log_section_header('Choosing k-mer range for assembly')

    # If the k-mer range file already exists, we use its values and proceed.
    kmer_range_filename = os.path.join(spades_dir, 'kmer_range')
    if os.path.isfile(kmer_range_filename):
        with open(kmer_range_filename, 'rt') as kmer_range_file:
            kmer_range = kmer_range_file.readline().strip().split(', ')
        if kmer_range:
            try:
                kmer_range = [int(x) for x in kmer_range]
                log.log('K-mer range already exists:')
                log.log('  ' + kmer_range_filename)
                log.log('\nWill use this existing range:')
                log.log('  ' + ', '.join([str(x) for x in kmer_range]))
                return kmer_range
            except ValueError:
                pass

    # If the code got here, then the k-mer range doesn't already exist and we'll create one by
    # examining the read lengths.
    read_lengths = get_read_lengths(reads_1_filename) + get_read_lengths(reads_2_filename)
    read_lengths = sorted(read_lengths)
    median_read_length = read_lengths[len(read_lengths) // 2 - 1]
    max_kmer = round_to_nearest_odd(max_kmer_frac * median_read_length)
    if max_kmer > 127:
        max_kmer = 127
    starting_kmer = round_to_nearest_odd(min_kmer_frac * max_kmer / max_kmer_frac)
    if starting_kmer < 11:
        starting_kmer = 11

    # Create the k-mer range from a non-linear function that spaces out the early k-mers more and
    # makes the later k-mers (which are most likely to be the good, used ones) closer together.
    kmer_range = []
    for x in [x / (kmer_count - 1) for x in range(kmer_count)]:
        kmer_range.append((max_kmer - starting_kmer) * (2 - 2 / (x + 1)) + starting_kmer)
    kmer_range = sorted(list(set([round_to_nearest_odd(x) for x in kmer_range])))
    kmer_range_str = ', '.join([str(x) for x in kmer_range])

    log.log('Median read length: ' + str(median_read_length))
    log.log('K-mer range: ' + kmer_range_str)

    kmer_range_file = open(kmer_range_filename, 'w')
    kmer_range_file.write(kmer_range_str)
    kmer_range_file.close()
    return kmer_range


def get_read_lengths(reads_filename):
    """
    Returns a list of the read lengths for the given read file.
    """
    if get_compression_type(reads_filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    reads = open_func(reads_filename, 'rb')
    read_lengths = []
    i = 0
    for line in reads:
        if i % 4 == 1:
            read_lengths.append(len(line.strip()))
        i += 1
    reads.close()
    return read_lengths


def count_segments_in_spades_fastg(fastg_file):
    seq_count = 0
    with open(fastg_file, 'rt') as fastg:
        for line in fastg:
            if line.startswith('>'):
                seq_count += 1
    return seq_count // 2
