"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module contains functionality related to miniasm, which Unicycler uses to build an assembly
using both Illumina contigs and long reads.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import os
import shutil
import statistics
import subprocess
import sys
from collections import defaultdict
from .misc import green, red, line_iterator, load_fasta, reverse_complement, print_table, \
    get_right_arrow
from .minimap_alignment import align_long_reads_to_assembly_graph, build_start_end_overlap_sets, \
    load_minimap_alignments_basic
from .string_graph import StringGraph, get_unsigned_seg_name
from . import log
from . import settings

try:
    from .cpp_wrappers import minimap_align_reads, miniasm_assembly, start_seq_alignment, \
        end_seq_alignment
except AttributeError as att_err:
    sys.exit('Error when importing C++ library: ' + str(att_err) + '\n'
             'Have you successfully build the library file using make?')


class MiniasmFailure(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


def build_miniasm_bridges(graph, out_dir, keep, threads, read_dict, long_read_filename,
                          scoring_scheme, racon_path):
    """
    EXTRACT READS USEFUL FOR LONG READ ASSEMBLY.
    * Take all single copy contigs over a certain length and get reads which overlap two or more.
      * While I'm at it, I should throw out reads which look like chimeras based on incompatible
        mapping.
    * Create a file of "long reads" which contains:
      * real long reads as found (and possibly split) by the above step
      * single copy contigs in FASTQ form (with a high quality, 'I' or something)

    """
    log.log_section_header('Assemble contigs and long reads with miniasm and Racon')
    log.log_explanation('Unicycler uses miniasm to construct a string graph '
                        'assembly using both the short read contigs and the long reads. If this '
                        'produces an assembly, Unicycler will extract bridges between '
                        'contigs, improve them with Racon and use them to simplify the assembly '
                        'graph. This method requires decent coverage of long reads and therefore '
                        'may not be fruitful if long reads are sparse. However, it does not '
                        'rely on the short read assembly graph having good connectivity and is '
                        'able to bridge an assembly graph even when it contains many dead ends.',
                        extra_empty_lines_after=0)
    log.log_explanation('Unicycler uses two types of "reads" as assembly input: '
                        'sufficiently long single-copy short read contigs and actual long reads '
                        'which overlap two or more of these contigs. It then assembles them with '
                        'a modified version of miniasm which gives precedence to the contigs over '
                        'the real long reads.', extra_empty_lines_after=0)
    log.log_explanation('Miniasm removes sequences which are contained in another sequence, and '
                        'this can result in short read contigs being lost in the string graph '
                        '(particularly common if the long reads are very long). If this happens, '
                        'the reads will be split where they contain a contig and Unicycler '
                        'will repeat the miniasm assembly until no contigs are lost.')

    miniasm_dir = os.path.join(out_dir, 'miniasm_assembly')
    if not os.path.exists(miniasm_dir):
        os.makedirs(miniasm_dir)

    # Align all the long reads to the graph and get the ones which overlap single-copy contigs
    # (and will therefore be useful for assembly).
    minimap_alignments = align_long_reads_to_assembly_graph(graph, long_read_filename,
                                                            miniasm_dir, threads)
    start_overlap_reads, end_overlap_reads = build_start_end_overlap_sets(minimap_alignments)
    assembly_read_names = get_miniasm_assembly_reads(minimap_alignments)

    assembly_reads_filename = os.path.join(miniasm_dir, '01_assembly_reads.fastq')
    mappings_filename = os.path.join(miniasm_dir, '02_mappings.paf')
    before_transitive_reduction_filename = os.path.join(miniasm_dir, '03_raw_string_graph.gfa')
    string_graph_filename = os.path.join(miniasm_dir, '10_final_string_graph.gfa')
    miniasm_output_filename = os.path.join(miniasm_dir, 'miniasm.out')

    # The miniasm assembly may be done multiple times. It reports whether a contig was contained in
    # a read, and if so, we will split the read and try again.
    read_break_points = defaultdict(set)
    mean_read_quals = {}
    for _ in range(20):  # limit to 20 iterations to prevent unforeseen infinite loop
        mean_read_quals = save_assembly_reads_to_file(assembly_reads_filename, assembly_read_names,
                                                      read_dict, graph, read_break_points)

        # Do an all-vs-all alignment of the assembly FASTQ, for miniasm input. Contig-contig
        # alignments are excluded (because single-copy contigs, by definition, should not overlap
        # each other).
        minimap_alignments_str = minimap_align_reads(assembly_reads_filename, assembly_reads_filename,
                                                     threads, 0, 'read vs read')
        with open(mappings_filename, 'wt') as mappings:
            for minimap_alignment_str in line_iterator(minimap_alignments_str):
                if minimap_alignment_str.count('CONTIG_') < 2:
                    mappings.write(minimap_alignment_str)
                    mappings.write('\n')

        # Now actually do the miniasm assembly, which will create a GFA file of the string graph.
        log.log('Assembling reads with miniasm... ', end='')
        min_depth = 3
        miniasm_assembly(assembly_reads_filename, mappings_filename, miniasm_dir, min_depth)

        contained_contig_count = 0
        with open(miniasm_output_filename, 'rt') as miniasm_out:
            for line in miniasm_out:
                line = line.strip()
                if line.startswith('CONTAINED CONTIG'):
                    line_parts = line.split('\t')
                    read_name = line_parts[2]
                    read_offset = 0
                    if '_range_' in read_name:
                        name_parts = read_name.split('_range_')
                        assert len(name_parts) == 2
                        read_name = name_parts[0]
                        read_offset = int(name_parts[1].split('-')[0])
                    read_start, read_end = int(line_parts[3]), int(line_parts[4])
                    read_break_points[read_name].add(read_offset + ((read_start + read_end) // 2))
                    contained_contig_count += 1

        # If the assembly finished without any contained contigs, then we're good to continue!
        if contained_contig_count == 0:
            break
        else:
            log.log(red('contigs lost in miniasm assembly'))
            log.log('Breaking reads and trying again\n')

    if not (os.path.isfile(string_graph_filename) and
            os.path.isfile(before_transitive_reduction_filename)):
        log.log(red('failed'))
        raise MiniasmFailure('miniasm failed to generate a string graph')
    string_graph = StringGraph(string_graph_filename, mean_read_quals)
    before_transitive_reduction = StringGraph(before_transitive_reduction_filename, mean_read_quals)

    log.log(green('success'))
    log.log('  ' + str(len(string_graph.segments)) + ' segments, ' +
            str(len(string_graph.links) // 2) + ' links', verbosity=2)

    string_graph.remove_branching_paths()
    if keep >= 3:
        string_graph.save_to_gfa(os.path.join(miniasm_dir, '12_branching_paths_removed.gfa'))

    string_graph.simplify_bridges(before_transitive_reduction)
    if keep >= 3:
        string_graph.save_to_gfa(os.path.join(miniasm_dir, '13_simplified_bridges.gfa'))

    string_graph.remove_overlaps(before_transitive_reduction, scoring_scheme)
    if keep >= 3:
        string_graph.save_to_gfa(os.path.join(miniasm_dir, '14_overlaps_removed.gfa'))

    string_graph.merge_reads()
    if keep >= 3:
        string_graph.save_to_gfa(os.path.join(miniasm_dir, '15_reads_merged.gfa'))

    # Some single copy contigs might be isolated from the main part of the graph (due to contained
    # read filtered or some graph simplification step, like bubble popping). We now need to place
    # them back in by aligning to the non-contig graph segments.
    string_graph.place_isolated_contigs(miniasm_dir, threads)
    if keep >= 3:
        string_graph.save_to_gfa(os.path.join(miniasm_dir, '18_contigs_placed.gfa'))

    # TO DO: I can probably remove this later, for efficiency. It's just a sanity check that none
    # of the graph manipulations screwed up the sequence ranges.
    string_graph.check_segment_names_and_ranges(read_dict, graph)

    # Polish each bridge sequence using Racon.
    polish_bridges(miniasm_dir, string_graph, read_dict, start_overlap_reads, end_overlap_reads,
                   racon_path, scoring_scheme, threads)
    string_graph.save_to_gfa(os.path.join(miniasm_dir, '19_racon_polish.gfa'))


    # TRY TO PLACE SMALLER SINGLE-COPY CONTIGS
    # * Use essentially the same process as place_isolated_contigs
    # *




    # LOOK FOR EACH BRIDGE SEQUENCE IN THE GRAPH.
    # * Goal 1: if we can find a short read version of the bridge, we should use that because it
    #   will probably be more accurate.
    # * Goal 2: using a graph path will let us 'use up' the segments, which helps with clean-up.
    # * In order to replace a miniasm assembly bridge sequence with a graph path sequence, the
    #   match has to be very strong! High identity over all sequence windows.
    # * Can use my existing path finding code, but tweak the settings to make them faster. This is
    #   because failing to find an existing path isn't too terrible, as we already have the miniasm
    #   sequence.

    # DO SOME BASIC GRAPH CLEAN-UP AND MERGE ALL POSSIBLE SEGMENTS.
    # * Clean up will be a bit tougher as we may have missed used sequence.

    # RE-RUN COPY NUMBER DETERMINATION.

    if keep < 3:
        shutil.rmtree(miniasm_dir)


def get_miniasm_assembly_reads(minimap_alignments):
    """
    Returns a list of read names which overlap at least one single copy graph segment.
    """
    miniasm_assembly_reads = []
    for read_name, alignments in minimap_alignments.items():
        if any(a.overlaps_reference() for a in alignments):
            miniasm_assembly_reads.append(read_name)
    return sorted(miniasm_assembly_reads)


def save_assembly_reads_to_file(read_filename, read_names, read_dict, graph, read_break_points):
    qual = chr(settings.CONTIG_READ_QSCORE + 33)
    log.log('Saving to ' + read_filename + ':')
    mean_read_quals = {}

    with open(read_filename, 'wt') as fastq:
        # First save the Illumina contigs as 'reads'. They are given a constant high qscore to
        # reflect our confidence in them.
        seg_count = 0
        for seg in sorted(graph.segments.values(), key=lambda x: x.number):
            if segment_suitable_for_miniasm_assembly(graph, seg):
                fastq.write('@CONTIG_')
                fastq.write(str(seg.number))
                fastq.write('\n')
                fastq.write(seg.forward_sequence)
                fastq.write('\n+\n')
                fastq.write(qual * seg.get_length())
                fastq.write('\n')
                seg_count += 1
        log.log('  ' + str(seg_count) + ' single copy contigs ' +
                str(settings.MIN_SEGMENT_LENGTH_FOR_MINIASM_BRIDGING) + ' bp or longer')

        # Now save the actual long reads.
        piece_count = 0
        for read_name in read_names:
            read = read_dict[read_name]
            if read_name not in read_break_points:
                breaks = [0, len(read.sequence)]
            else:
                breaks = [0] + sorted(read_break_points[read_name]) + [len(read.sequence)]
            read_ranges = [(breaks[i-1], breaks[i]) for i in range(1, len(breaks))]

            for i, read_range in enumerate(read_ranges):
                s, e = read_range
                seq = read.sequence[s:e]
                if len(seq) < 100:
                    continue
                quals = read.qualities[s:e]
                if len(read_ranges) == 1:
                    read_range_name = read_name
                else:
                    read_range_name = read_name + '_range_' + str(s) + '-' + str(e)
                fastq.write('@')
                fastq.write(read_range_name)
                fastq.write('\n')
                fastq.write(seq)
                fastq.write('\n+\n')
                fastq.write(quals)
                fastq.write('\n')
                mean_read_quals[read_range_name] = statistics.mean(ord(x)-33 for x in quals)
                piece_count += 1

        break_string = ''
        if piece_count > len(read_names):
            break_string = ' broken into ' + str(piece_count) + ' pieces'
        log.log('  ' + str(len(read_names)) + ' overlapping long reads (out of ' +
                str(len(read_dict)) + ' total long reads)' + break_string)
    log.log('')
    return mean_read_quals


def segment_suitable_for_miniasm_assembly(graph, segment):
    """
    Returns True if the segment is:
      1) single copy
      2) long enough
      3) not already circular and complete
    """
    if graph.get_copy_number(segment) != 1:
        return False
    if segment.get_length() < settings.MIN_SEGMENT_LENGTH_FOR_MINIASM_BRIDGING:
        return False
    return not graph.is_component_complete([segment.number])


def polish_bridges(miniasm_dir, string_graph, read_dict, start_overlap_reads, end_overlap_reads,
                   racon_path, scoring_scheme, threads):
    log.log_explanation('Unicycler now uses Racon to polish each string graph segment which '
                        'connects two contigs. The resulting consensus sequences will be more '
                        'accurate than the starting sequences (which are made from single long '
                        'reads).')
    ra = get_right_arrow()
    col_widths = [22, len(ra), 22, 6, 5, 10, 10]
    racon_table_header = ['Start', '', 'End', 'Racon rounds', 'Read depth', 'Pre-Racon length',
                          'Post-Racon length']
    print_table([racon_table_header], fixed_col_widths=col_widths, left_align_header=False,
                alignments='RLLRRRR', indent=0, col_separation=1)

    bridging_paths = string_graph.get_bridging_paths()
    for bridging_path in bridging_paths:
        assert len(bridging_path) == 3
        pre_racon_length = len(string_graph.seq_from_signed_seg_name(bridging_path[1]))
        rounds, depth = polish_bridge(bridging_path, miniasm_dir, string_graph, read_dict,
                                      start_overlap_reads, end_overlap_reads, racon_path,
                                      scoring_scheme, threads)
        post_racon_length = len(string_graph.seq_from_signed_seg_name(bridging_path[1]))
        racon_table_row = [bridging_path[0], ra, bridging_path[2], str(rounds), '%.1f' % depth,
                           str(pre_racon_length), str(post_racon_length)]
        print_table([racon_table_row], fixed_col_widths=col_widths, left_align_header=False,
                    alignments='RLLRRRR', indent=0, header_format='normal', col_separation=1)


def polish_bridge(bridging_path, miniasm_dir, string_graph, read_dict, start_overlap_reads,
                  end_overlap_reads, racon_path, scoring_scheme, threads):
    first = string_graph.segments[get_unsigned_seg_name(bridging_path[0])]
    middle = string_graph.segments[get_unsigned_seg_name(bridging_path[1])]
    last = string_graph.segments[get_unsigned_seg_name(bridging_path[2])]
    first_sign = bridging_path[0][-1]
    middle_sign = bridging_path[1][-1]
    last_sign = bridging_path[2][-1]
    first_num = int(first.short_name.split('CONTIG_')[-1]) * (-1 if first_sign == '-' else 1)
    last_num = int(last.short_name.split('CONTIG_')[-1]) * (-1 if last_sign == '-' else 1)

    bridge_name = 'bridge_' + (first.short_name + '_to_' + last.short_name).replace('CONTIG_', '')
    bridge_dir = os.path.join(miniasm_dir, bridge_name)
    if not os.path.isdir(bridge_dir):
        os.makedirs(bridge_dir)

    # Save the bridge to file, with a bit of contig sequence on each end.
    margin = settings.RACON_POLISH_MARGIN
    bridge_start_seq = string_graph.seq_from_signed_seg_name(bridging_path[0])[-margin:]
    bridge_middle_seq = string_graph.seq_from_signed_seg_name(bridging_path[1])
    bridge_end_seq = string_graph.seq_from_signed_seg_name(bridging_path[2])[:margin]
    full_bridge_seq = bridge_start_seq + bridge_middle_seq + bridge_end_seq
    unpolished_seq = os.path.join(bridge_dir, 'unpolished.fasta')
    with open(unpolished_seq, 'wt') as fasta:
        fasta.write('>bridge\n')
        fasta.write(full_bridge_seq)
        fasta.write('\n')

    # Save the reads to file which we'll use to polish.
    read_names = sorted(end_overlap_reads[first_num] | start_overlap_reads[last_num])
    if not read_names:
        return 0
    polish_reads = os.path.join(bridge_dir, 'reads.fastq')
    with open(polish_reads, 'wt') as fastq:
        for name in read_names:
            fastq.write(read_dict[name].get_fastq())

    current_fasta = unpolished_seq
    last_seq = full_bridge_seq
    previous_seqs = [last_seq]

    polish_round_count = 0
    depth = 0.0
    for _ in range(settings.RACON_POLISH_LOOP_COUNT):
        round_num_str = '%02d' % (polish_round_count + 1)

        # Create the alignments for Racon using miniasm.
        mappings_filename = os.path.join(bridge_dir, 'alignments_' + round_num_str + '.paf')
        minimap_alignments_str = minimap_align_reads(current_fasta, polish_reads, 1, 0)
        with open(mappings_filename, 'wt') as mappings:
            mappings.write(minimap_alignments_str)
        depth = sum(a.fraction_ref_aligned()
                    for a in load_minimap_alignments_basic(minimap_alignments_str))

        # Run Racon. It crashes sometimes, so repeat until its return code is 0.
        racon_fasta = os.path.join(bridge_dir, 'consensus_' + round_num_str + '.fasta')
        racon_log = os.path.join(bridge_dir, 'consensus_' + round_num_str + '.log')
        command = [racon_path, '--verbose', '9', '--threads', str(threads), '--bq', '-1',
                   polish_reads, mappings_filename, current_fasta, racon_fasta]
        return_code = 1
        for t in range(100):  # Only try a fixed number of times, to prevent an infinite loop.
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = process.communicate()
            with open(racon_log, 'wb') as log_file:
                log_file.write(out)
                log_file.write(err)
            return_code = process.returncode
            if return_code == 0:
                break
            if return_code != 0:
                if os.path.isfile(racon_fasta):
                    os.remove(racon_fasta)
                if os.path.isfile(racon_log):
                    os.remove(racon_log)
        if return_code != 0 or not os.path.isfile(racon_fasta):
            break
        polished_seq = load_fasta(racon_fasta)[0][1]
        polish_round_count += 1

        # If the sequence hasn't changed or is equal to any previous state (sometimes Racon
        # oscillates between a few states), there's no need to do another round of polishing.
        if any(polished_seq == x for x in previous_seqs):
            break

        current_fasta = racon_fasta
        last_seq = polished_seq
        previous_seqs.append(last_seq)

    # Use some alignments to figure out where the contigs are in the polished sequence so we can
    # snip out the bridge part of the sequence and save it back into the graph segment.
    middle_start = start_seq_alignment(bridge_start_seq, last_seq, scoring_scheme)
    middle_end = end_seq_alignment(bridge_end_seq, last_seq, scoring_scheme)
    new_bridge_seq = last_seq[middle_start:middle_end]
    new_rev_bridge_seq = reverse_complement(new_bridge_seq)
    if middle_sign == '-':
        new_bridge_seq, new_rev_bridge_seq = new_rev_bridge_seq, new_bridge_seq
    middle.forward_sequence = new_bridge_seq
    middle.reverse_sequence = new_rev_bridge_seq

    return polish_round_count, depth
