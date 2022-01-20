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
from .misc import round_to_nearest_odd, get_compression_type, int_to_str, quit_with_error, \
    bold, dim, print_table, get_left_arrow, float_to_str
from .assembly_graph import AssemblyGraph
from . import log


class BadFastq(Exception):
    pass


def get_best_spades_graph(short1, short2, short_unpaired, out_dir, read_depth_filter, verbosity,
                          spades_path, threads, keep, kmer_count, min_k_frac, max_k_frac, kmers,
                          expected_linear_seqs, largest_component, spades_graph_prefix,
                          spades_options):
    """
    This function tries a SPAdes assembly at different k-mers and returns the best one.
    """
    spades_dir = os.path.join(out_dir, 'spades_assembly')
    if not os.path.exists(spades_dir):
        os.makedirs(spades_dir)

    threads = min(threads, 32)  # SPAdes can possibly crash if given too many threads.
    check_fastqs(short1, short2, short_unpaired)
    reads = (short1, short2, short_unpaired)

    kmer_range = get_kmer_range(kmers, short1, short2, short_unpaired, spades_dir, kmer_count,
                                min_k_frac, max_k_frac, spades_path)

    log.log_section_header('SPAdes assemblies')
    log.log_explanation('Unicycler now uses SPAdes to assemble the short reads. It scores the '
                        'assembly graph for each k-mer using the number of contigs (fewer is '
                        'better) and the number of dead ends (fewer is better). The score '
                        'function is 1/(c*(d+2)), where c is the contig count and d is the '
                        'dead end count.')

    if verbosity > 1:
        spades_results_table = [['K-mer', 'Contigs', 'Links', 'Total length', 'N50',
                                 'Longest contig', 'Dead ends', 'Score']]
    else:
        spades_results_table = [['K-mer', 'Contigs', 'Dead ends', 'Score']]

    graph_files, insert_size_mean, insert_size_deviation = \
        run_spades_all_kmers(reads, spades_dir, kmer_range, threads, spades_path,
                             spades_graph_prefix, spades_options)

    existing_graph_files = [x for x in graph_files if x is not None]
    if not existing_graph_files:
        quit_with_error('SPAdes failed to produce assemblies. '
                        'See spades_assembly/spades.log for more info.')
    median_segment_count = statistics.median(count_segments_in_gfa(x)
                                             for x in existing_graph_files)

    best_score, best_kmer, best_graph_filename = 0.0, 0, ''
    for graph_file, kmer in zip(graph_files, kmer_range):
        table_line = [int_to_str(kmer)]

        if graph_file is None:
            table_line += [''] * (7 if verbosity > 1 else 2)
            table_line.append('failed')
            spades_results_table.append(table_line)
            continue

        assembly_graph = AssemblyGraph(graph_file, kmer, insert_size_mean=insert_size_mean,
                                       insert_size_deviation=insert_size_deviation)

        # If this graph has way too many segments, then we will just skip it because very complex
        # graphs take forever to clean up.
        # TO DO: I can remove this awkward hack if I make the graph cleaning more efficient.
        if len(assembly_graph.segments) > 4 * median_segment_count:
            table_line += [''] * (6 if verbosity > 1 else 2)
            table_line.append('too complex')
            spades_results_table.append(table_line)
            continue

        log.log('\nCleaning k{} graph'.format(kmer), 2)
        assembly_graph.clean(read_depth_filter, largest_component)
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
            score = 1.0 / (segment_count * (adjusted_dead_ends + 2))

        # Prepare the table line for this k-mer graph.
        table_line += [int_to_str(segment_count)]
        if verbosity > 1:
            n50, shortest, _, median, _, longest = assembly_graph.get_contig_stats()
            table_line += [int_to_str(assembly_graph.get_total_link_count()),
                           int_to_str(assembly_graph.get_total_length()),
                           int_to_str(n50), int_to_str(longest)]
        table_line += [int_to_str(dead_ends), '{:.2e}'.format(score)]
        spades_results_table.append(table_line)

        if score > best_score:
            best_kmer, best_score, best_graph_filename = kmer, score, graph_file

    log.log('', 2)

    if not best_kmer:
        quit_with_error('none of the SPAdes graphs were suitable for scaffolding in Unicycler')

    assembly_graph = AssemblyGraph(best_graph_filename, best_kmer,
                                   insert_size_mean=insert_size_mean,
                                   insert_size_deviation=insert_size_deviation)
    removed_count, removed_length = assembly_graph.clean(read_depth_filter, largest_component)
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

    # Report on the results of the read depth filter (can help with identifying levels of
    # contamination).
    log.log('\nRead depth filter: removed {} contigs totalling {} bp'.format(removed_count,
                                                                             removed_length))
    # Clean up.
    if keep < 1:
        for g in graph_files:
            if os.path.isfile(g):
                log.log('Deleting ' + g)
                os.remove(g)
    if keep < 3 and os.path.isdir(spades_dir):
        log.log('Deleting ' + spades_dir + '/')
        shutil.rmtree(spades_dir, ignore_errors=True)

    return assembly_graph


def run_spades_all_kmers(read_files, spades_dir, kmers, threads, spades_path, spades_graph_prefix,
                         spades_options):
    """
    SPAdes is run with all k-mers up to the top one. For example:
      * round 1: 25
      * round 2: 25,43
      * round 3: 25,43,59
      * round 4: 25,43,59,73

    This is because it only saves the necessary graph information for the final k-mer. The first
    round is a normal SPAdes run, and subsequent rounds use the --restart-from option.
    """
    short1, short2, unpaired = read_files[0], read_files[1], read_files[2]
    using_paired_reads = short1 is not None and short2 is not None and \
        os.path.isfile(short1) and os.path.isfile(short2)
    using_unpaired_reads = unpaired is not None and os.path.isfile(unpaired)

    graph_files, insert_size_means, insert_size_deviations = [], [], []
    for i in range(len(kmers)):
        biggest_kmer = kmers[i]
        command = build_spades_command(spades_path, spades_dir, threads, kmers, i, short1, short2,
                                       unpaired, using_paired_reads, using_unpaired_reads,
                                       spades_options)
        log.log(' '.join(command))
        graph_file, insert_size_mean, insert_size_deviation = \
            run_spades_one_kmer(command, spades_dir, biggest_kmer)

        copy_path = spades_graph_prefix + '_K' + '{:03d}'.format(biggest_kmer) + '.gfa'
        shutil.copy(graph_file, copy_path)
        graph_files.append(copy_path)
        insert_size_means.append(insert_size_mean)
        insert_size_deviations.append(insert_size_deviation)
    log.log('')

    insert_size_means = [x for x in insert_size_means if x is not None]
    insert_size_deviations = [x for x in insert_size_deviations if x is not None]

    if insert_size_means and insert_size_deviations:
        insert_size_mean = statistics.median(insert_size_means)
        insert_size_deviation = statistics.median(insert_size_deviations)
    else:
        # If we couldn't get the insert size from the SPAdes output (e.g. it was an
        # unpaired-reads-only assembly), we'll use the read length instead.
        read_lengths = get_read_lengths(short1) + get_read_lengths(short2) + \
                       get_read_lengths(unpaired)
        insert_size_mean = statistics.mean(read_lengths)
        insert_size_deviation = max(statistics.stdev(read_lengths), 1.0)

    log.log('', 2)
    log.log('Insert size mean: ' + float_to_str(insert_size_mean, 1) + ' bp', 2)
    log.log('Insert size stdev: ' + float_to_str(insert_size_deviation, 1) + ' bp', 2)
    log.log('', 2)

    return graph_files, insert_size_mean, insert_size_deviation


def build_spades_command(spades_path, spades_dir, threads, kmers, i, short1, short2, unpaired,
                         using_paired_reads, using_unpaired_reads, spades_options):
    kmer_string = ','.join([str(x) for x in kmers[:i+1]])

    command = [spades_path, '-o', spades_dir, '-k', kmer_string, '--threads', str(threads)]
    if i == 0:  # first k-mer
        command += ['--isolate']
        if using_paired_reads:
            command += ['-1', short1, '-2', short2]
        if using_unpaired_reads:
            command += ['-s', unpaired]
    else:  # subsequent k-mer
        previous_k = kmers[i - 1]
        command += ['--restart-from', f'k{previous_k}']
    if spades_options:
        command += spades_options.split()
    if not spades_options or '-m' not in spades_options.split():
        command += ['-m', '1024']
    return command


def run_spades_one_kmer(command, spades_dir, biggest_kmer):
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

        try:
            insert_size_mean = float(spades_output.split('Insert size = ')[-1].split(',')[0])
            insert_size_deviation = float(spades_output.split('deviation = ')[-1].split(',')[0])
        except ValueError:
            pass

    if process.returncode != 0:
        spades_error = process.stderr.read().strip().decode()
        quit_with_error('SPAdes encountered an error:\n' + spades_error)

    graph_file = os.path.join(spades_dir, 'K' + str(biggest_kmer),
                              'assembly_graph_with_scaffolds.gfa')
    if not os.path.isfile(graph_file):
        quit_with_error('SPAdes failed to produce an assembly graph: ' + str(graph_file))

    return graph_file, insert_size_mean, insert_size_deviation


def check_fastqs(short1, short2, short_unpaired):
    using_paired_reads = bool(short1) and bool(short2)
    using_unpaired_reads = bool(short_unpaired)
    if using_paired_reads:
        count_1, count_2 = 0, 0
        try:
            count_1 = get_read_count(short1)
        except BadFastq:
            quit_with_error('this read file is not a properly formatted FASTQ: ' + short1)
        try:
            count_2 = get_read_count(short2)
        except BadFastq:
            quit_with_error('this read file is not a properly formatted FASTQ: ' + short2)
        if count_1 != count_2:
            quit_with_error('the paired read input files have an unequal number of reads')
    if using_unpaired_reads:
        try:
            get_read_count(short_unpaired)
        except BadFastq:
            quit_with_error('this read file is not properly formatted as FASTQ: ' + short_unpaired)


def get_max_spades_kmer(spades_path):
    """
    SPAdes usually has a maximum k-mer size of 127, but this can be changed when compiling SPAdes,
    so this function checks the help text to see what it is.
    https://github.com/ablab/spades/issues/40
    """
    try:
        process = subprocess.Popen([spades_path, '--help'], stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        out, err = process.communicate()
        all_output = out.decode() + err.decode()
        all_output = all_output.replace('\n', ' ')
        all_output = ' '.join(all_output.split())
        max_kmer = all_output.split('must be odd and less than ')[1].split(')')[0]
        return int(max_kmer) - 1
    except (IndexError, ValueError):
        return 127


def get_kmer_range(given_kmers, reads_1_filename, reads_2_filename, unpaired_reads_filename,
                   spades_dir, kmer_count, min_kmer_frac, max_kmer_frac, spades_path):
    """
    Uses the read lengths to determine the k-mer range to be used in the SPAdes assembly.
    """
    # If the user specified which k-mers to use, we can skip automatic k-mer selection.
    if given_kmers is not None:
        return given_kmers

    log.log_section_header('Choosing k-mer range for assembly')
    log.log_explanation('Unicycler chooses a k-mer range for SPAdes based on the length of the '
                        'input reads. It uses a wide range of many k-mer sizes to maximise the '
                        'chance of finding an ideal assembly.')

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

    max_spades_kmer = get_max_spades_kmer(spades_path)
    log.log('SPAdes maximum k-mer: {}'.format(max_spades_kmer))
    if max_spades_kmer != 127:
        log.log('    (unusual value, probably indicates custom SPAdes compilation)')

    # If the code got here, then the k-mer range doesn't already exist and we'll create one by
    # examining the read lengths.
    read_lengths = get_read_lengths(reads_1_filename) + get_read_lengths(reads_2_filename) + \
        get_read_lengths(unpaired_reads_filename)
    read_lengths = sorted(read_lengths)
    median_read_length = read_lengths[len(read_lengths) // 2 - 1]
    max_kmer = round_to_nearest_odd(max_kmer_frac * median_read_length)
    if max_kmer > max_spades_kmer:
        max_kmer = max_spades_kmer
    starting_kmer = round_to_nearest_odd(min_kmer_frac * max_kmer / max_kmer_frac)
    if starting_kmer < 11:
        starting_kmer = 11

    if kmer_count == 1:
        kmer_range = [max_kmer]
    elif kmer_count == 2:
        kmer_range = [starting_kmer, max_kmer]
    else:
        # Create the k-mer range from a non-linear function that spaces out the early k-mers more
        # and makes the later k-mers (which are most likely to be the good, used ones) closer
        # together.
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
    if reads_filename is None:
        return []
    if get_compression_type(reads_filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    with open_func(reads_filename, 'rb') as reads:
        read_lengths = []
        i = 0
        for line in reads:
            if i % 4 == 1:
                read_lengths.append(len(line.strip()))
            i += 1
    return read_lengths


def get_read_count(reads_filename):
    """
    Returns the number of reads in the given file.
    """
    if reads_filename is None:
        return 0
    if get_compression_type(reads_filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    with open_func(reads_filename, 'rb') as reads:
        read_count = 0
        i = 0
        for line in reads:
            if i % 4 == 0:
                try:
                    assert line.startswith(b'@')
                except AssertionError:
                    raise BadFastq
                read_count += 1
            i += 1
    return read_count


def count_segments_in_gfa(fastg_file):
    seq_count = 0
    with open(fastg_file, 'rt') as fastg:
        for line in fastg:
            if line.startswith('S\t'):
                seq_count += 1
    return seq_count
