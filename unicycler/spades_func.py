"""
Functions relating to SPAdes assembly steps.

Author: Ryan Wick
email: rrwick@gmail.com
"""

import os
import subprocess
import gzip
import shutil
import statistics
from .misc import print_section_header, round_to_nearest_odd, get_compression_type, int_to_str, \
    quit_with_error, strip_read_extensions, bold, dim, print_table, print_v, get_timestamp, \
    get_left_arrow
from .assembly_graph import AssemblyGraph


def get_best_spades_graph(short1, short2, short_unpaired, out_dir, read_depth_filter, verbosity,
                          spades_path, threads, keep_temp, kmer_count, min_kmer_frac, max_kmer_frac,
                          no_spades_correct, expected_linear_seqs):
    """
    This function tries a SPAdes assembly at different kmers and returns the best.
    'The best' is defined as the smallest dead-end count after low-depth filtering.  If multiple
    graphs have the same dead-end count (e.g. zero!) then the highest kmer is used.
    """
    spades_dir = os.path.join(out_dir, 'spades_assembly_temp')
    if not os.path.exists(spades_dir):
        os.makedirs(spades_dir)

    if no_spades_correct:
        reads = (short1, short2, short_unpaired)
    else:
        reads = spades_read_correction(short1, short2, short_unpaired, spades_dir, verbosity,
                                       threads, spades_path)
    kmer_range = get_kmer_range(short1, short2, spades_dir, verbosity, kmer_count, min_kmer_frac,
                                max_kmer_frac)
    assem_dir = os.path.join(spades_dir, 'assembly')
    print_section_header('Conducting SPAdes assemblies', verbosity)

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
        spades_assembly(reads, assem_dir, kmer_range, verbosity, threads, spades_path)

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

        print_v(get_timestamp() + ' cleaning SPAdes graph (k-mer=' + str(kmer) + ')', verbosity, 3)
        assembly_graph.clean(read_depth_filter, verbosity)
        clean_graph_filename = os.path.join(spades_dir, 'k' + str(kmer) + '_assembly_graph.gfa')
        assembly_graph.save_to_gfa(os.path.join(spades_dir, clean_graph_filename), 0)

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

    if not best_kmer:
        quit_with_error('none of the SPAdes graphs were suitable for scaffolding in Unicycler')

    # If the best k-mer is the top k-mer, then SPAdes has already done the repeat resolution and
    # we can just grab it now. Easy! If the best k-mer was a different k-mer size, then we need
    # to run SPAdes again to get that repeat resolution.
    if best_kmer != kmer_range[-1]:
        new_kmer_range = [x for x in kmer_range if x <= best_kmer]
        graph_file, insert_size_mean, insert_size_deviation = \
            spades_assembly(reads, assem_dir, new_kmer_range, verbosity, threads, spades_path,
                            just_last=True)
        best_graph_filename = graph_file
    paths_file = os.path.join(assem_dir, 'contigs.paths')
    if os.path.isfile(paths_file):
        copied_paths_file = os.path.join(spades_dir,
                                         'k' + ('%03d' % best_kmer) + '_contigs.paths')
        shutil.copyfile(paths_file, copied_paths_file)

    # Now we can load and clean the graph again, this time giving it the SPAdes contig paths.
    assembly_graph = AssemblyGraph(best_graph_filename, best_kmer, paths_file=paths_file,
                                   insert_size_mean=insert_size_mean,
                                   insert_size_deviation=insert_size_deviation)
    assembly_graph.clean(read_depth_filter, verbosity)
    clean_graph_filename = os.path.join(spades_dir,
                                        'k' + str(best_kmer) + '_assembly_graph.gfa')
    assembly_graph.save_to_gfa(os.path.join(spades_dir, clean_graph_filename), 0)

    # Clean up SPAdes files.
    if keep_temp < 2 and os.path.isdir(assem_dir):
        shutil.rmtree(assem_dir)
    if keep_temp < 1 and os.path.isdir(spades_dir):
        shutil.rmtree(spades_dir)

    if best_score == 0.0:
        quit_with_error('none of the SPAdes assemblies produced assembled sequence')

    # Print the SPAdes result table, highlighting the best k-mer in green.
    if verbosity > 1:
        print_section_header('SPAdes assembly graph summary', verbosity)
    if verbosity > 0:
        best_kmer_row = [x[0] for x in spades_results_table].index(int_to_str(best_kmer))
        print_table(spades_results_table, alignments='RRRRRRRR', indent=0,
                    row_colour={best_kmer_row: 'green'},
                    row_extra_text={best_kmer_row: ' ' + get_left_arrow() + 'best'})

    return assembly_graph


def spades_read_correction(short1, short2, short_unpaired, spades_dir, verbosity, threads,
                           spades_path):
    """
    This runs SPAdes with the --only-error-correction option.
    """
    print_section_header('SPAdes read error correction', verbosity)

    # If the corrected reads already exist, then we just use them and proceed.
    corrected_1 = os.path.join(spades_dir, 'corrected_1.fastq.gz')
    corrected_2 = os.path.join(spades_dir, 'corrected_2.fastq.gz')
    corrected_u = os.path.join(spades_dir, 'corrected_u.fastq.gz')
    if os.path.isfile(corrected_1) and os.path.isfile(corrected_2):
        if verbosity > 0:
            print('Corrected reads already exist. Will use these reads instead of running SPAdes '
                  'error correction:')
            print('  ' + corrected_1)
            print('  ' + corrected_2)
            if os.path.isfile(corrected_u):
                print('  ' + corrected_u)
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
        if verbosity > 0 and 'Command line:' in spades_output:
            spades_output = ' '.join(spades_output.split()).replace('Command line: ', '')
            print('Command: ' + bold(spades_output), flush=True)
        elif spades_output and verbosity > 1:
            print(dim(spades_output), flush=True)

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

    if verbosity > 0:
        print()
        print('Corrected reads:')
        print('  ' + corrected_1)
        print('  ' + corrected_2)
        if corrected_u:
            print('  ' + corrected_u)

    return corrected_1, corrected_2, corrected_u


def spades_assembly(read_files, out_dir, kmers, verbosity, threads, spades_path, just_last=False):
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
        if spades_output and verbosity > 1:
            if spades_output.startswith('Command line:'):
                spades_output = ' '.join(spades_output.split()).replace('Command line: ', '')
                print('Command: ' + bold(spades_output))
            else:
                print(dim(spades_output), flush=True)
        if 'Insert size =' in spades_output and 'deviation = ' in spades_output:
            try:
                insert_size_mean = float(spades_output.split('Insert size = ')[-1].split(',')[0])
                insert_size_deviation = float(spades_output.split('deviation = ')[-1].split(',')[0])
            except ValueError:
                pass

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
                                                 'k' + ('%03d' % kmer) + '_assembly_graph.fastg')
                shutil.copyfile(graph_file, copied_graph_file)
                graph_files.append(copied_graph_file)
            else:
                graph_files.append(None)
        return graph_files, insert_size_mean, insert_size_deviation


def get_kmer_range(reads_1_filename, reads_2_filename, spades_dir, verbosity, kmer_count,
                   min_kmer_frac, max_kmer_frac):
    """
    Uses the read lengths to determine the k-mer range to be used in the SPAdes assembly.
    """
    print_section_header('Choosing k-mer range for assembly', verbosity)

    # If the k-mer range file already exists, we use its values and proceed.
    kmer_range_filename = os.path.join(spades_dir, 'kmer_range')
    if os.path.isfile(kmer_range_filename):
        with open(kmer_range_filename, 'rt') as kmer_range_file:
            kmer_range = kmer_range_file.readline().strip().split(', ')
        if kmer_range:
            try:
                kmer_range = [int(x) for x in kmer_range]
                if verbosity > 0:
                    print('K-mer range already exists:')
                    print('  ' + kmer_range_filename)
                    print('\nWill use this existing range:')
                    print('  ' + ', '.join([str(x) for x in kmer_range]))
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

    if verbosity > 0:
        print('Median read length: ' + str(median_read_length))
        print('K-mer range: ' + kmer_range_str)

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
