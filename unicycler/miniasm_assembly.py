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
from .misc import green, red, line_iterator, print_table, int_to_str, remove_dupes_preserve_order, \
    reverse_complement
from .minimap_alignment import align_long_reads_to_assembly_graph, \
    load_minimap_alignments_basic, load_minimap_alignments, combine_close_hits, \
    remove_conflicting_alignments
from .string_graph import StringGraph, merge_string_graph_segments_into_unitig_graph
from . import log
from . import settings

try:
    from .cpp_wrappers import minimap_align_reads, miniasm_assembly, start_seq_alignment, \
        end_seq_alignment
except AttributeError as att_err:
    sys.exit('Error when importing C++ library: ' + str(att_err) + '\n'
             'Have you successfully built the library file using make?')


class MiniasmFailure(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


def make_miniasm_string_graph(graph, out_dir, keep, threads, read_dict, long_read_filename,
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
    log.log_section_header('Assembling contigs and long reads with miniasm')
    log.log_explanation('Unicycler uses miniasm to construct a string graph '
                        'assembly using both the short read contigs and the long reads. It can '
                        'then use the resulting string graph to produce bridges between contigs. '
                        'This method requires decent coverage of long reads and therefore '
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
                        '(particularly if the long reads are very long). If this happens, '
                        'the reads will be split where they contain a contig and Unicycler '
                        'will repeat the miniasm assembly until no contigs are lost.')

    miniasm_dir = os.path.join(out_dir, 'miniasm_assembly')
    if not os.path.exists(miniasm_dir):
        os.makedirs(miniasm_dir)

    # Align all the long reads to the graph and get the ones which overlap single-copy contigs
    # (and will therefore be useful for assembly).
    minimap_alignments = align_long_reads_to_assembly_graph(graph, long_read_filename,
                                                            miniasm_dir, threads)
    assembly_read_names = get_miniasm_assembly_reads(minimap_alignments)

    assembly_reads_filename = os.path.join(miniasm_dir, '01_assembly_reads.fastq')
    mappings_filename = os.path.join(miniasm_dir, '02_mappings.paf')
    before_transitive_reduction_filename = os.path.join(miniasm_dir, '03_raw_string_graph.gfa')
    string_graph_filename = os.path.join(miniasm_dir, '10_final_string_graph.gfa')
    branching_paths_removed_filename = os.path.join(miniasm_dir, '11_branching_paths_removed.gfa')
    unitig_graph_filename = os.path.join(miniasm_dir, '12_unitig_graph.gfa')
    racon_polished_filename = os.path.join(miniasm_dir, '13_racon_polished.gfa')
    contigs_placed_filename = os.path.join(miniasm_dir, '14_contigs_placed.gfa')

    mean_read_quals = save_assembly_reads_to_file(assembly_reads_filename, assembly_read_names,
                                                  read_dict, graph)

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

    # TO DO: refine these overlaps? Perhaps using Unicycler-align? I suspect that the quality of a
    # miniasm assembly is highly dependent on the quality of the input overlaps.

    # Now actually do the miniasm assembly, which will create a GFA file of the string graph.
    log.log('Assembling reads with miniasm... ', end='')
    min_depth = 3
    miniasm_assembly(assembly_reads_filename, mappings_filename, miniasm_dir, min_depth)
    if not (os.path.isfile(string_graph_filename) and
            os.path.isfile(before_transitive_reduction_filename)):
        log.log(red('failed'))
        raise MiniasmFailure('miniasm failed to generate a string graph')
    string_graph = StringGraph(string_graph_filename, mean_read_quals)

    log.log(green('success'))
    log.log('  ' + str(len(string_graph.segments)) + ' segments, ' +
            str(len(string_graph.links) // 2) + ' links', verbosity=2)

    string_graph.remove_branching_paths()
    if keep >= 3:
        string_graph.save_to_gfa(branching_paths_removed_filename)

    unitig_graph = merge_string_graph_segments_into_unitig_graph(string_graph)
    if keep >= 3:
        unitig_graph.save_to_gfa(unitig_graph_filename)

    polish_unitigs_with_racon(unitig_graph, miniasm_dir, assembly_read_names, read_dict, graph,
                              racon_path, threads, scoring_scheme)
    if keep >= 3:
        unitig_graph.save_to_gfa(racon_polished_filename)

    place_contigs(miniasm_dir, graph, unitig_graph, threads)
    if keep >= 3:
        unitig_graph.save_to_gfa(contigs_placed_filename)














    quit()  # TEMP




























    if keep < 3:
        shutil.rmtree(miniasm_dir)
    return string_graph


def get_miniasm_assembly_reads(minimap_alignments):
    """
    Returns a list of read names which overlap at least one single copy graph segment.
    """
    miniasm_assembly_reads = []
    for read_name, alignments in minimap_alignments.items():
        if any(a.overlaps_reference() for a in alignments):
            miniasm_assembly_reads.append(read_name)
    return sorted(miniasm_assembly_reads)


def save_assembly_reads_to_file(read_filename, read_names, read_dict, graph, contig_copy_count=1):
    qual = chr(settings.CONTIG_READ_QSCORE + 33)
    log.log('Saving to ' + read_filename + ':')
    mean_read_quals = {}

    with open(read_filename, 'wt') as fastq:
        # First save the Illumina contigs as 'reads'. They are given a constant high qscore to
        # reflect our confidence in them.
        seg_count = 0
        for seg in sorted(graph.segments.values(), key=lambda x: x.number):
            if segment_suitable_for_miniasm_assembly(graph, seg):
                for i in range(contig_copy_count):
                    fastq_name = '@CONTIG_' + str(seg.number)
                    if contig_copy_count > 1:
                        fastq_name += '_' + str(i+1)
                    fastq.write(fastq_name)
                    fastq.write('\n')
                    if i % 2 == 0:  # evens
                        fastq.write(seg.forward_sequence)
                    else:  # odds
                        fastq.write(seg.reverse_sequence)
                    fastq.write('\n+\n')
                    fastq.write(qual * seg.get_length())
                    fastq.write('\n')
                seg_count += 1
        message = '  ' + str(seg_count) + ' single copy contigs ' + \
                  str(settings.MIN_SEGMENT_LENGTH_FOR_MINIASM_BRIDGING) + ' bp or longer'
        if contig_copy_count > 1:
            message += ' (duplicated ' + str(contig_copy_count) + ' times)'
        log.log(message)

        # Now save the actual long reads.
        for read_name in read_names:
            read = read_dict[read_name]
            seq = read.sequence
            if len(seq) < 100:
                continue
            quals = read.qualities
            fastq.write('@')
            fastq.write(read_name)
            fastq.write('\n')
            fastq.write(seq)
            fastq.write('\n+\n')
            fastq.write(quals)
            fastq.write('\n')
            mean_read_quals[read_name] = statistics.mean(ord(x)-33 for x in quals)

        log.log('  ' + str(len(read_names)) + ' overlapping long reads (out of ' +
                str(len(read_dict)) + ' total long reads)')
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


def polish_unitigs_with_racon(unitig_graph, miniasm_dir, assembly_read_names, read_dict, graph,
                              racon_path, threads, scoring_scheme):
    log.log_section_header('Polishing miniasm assembly with Racon')
    log.log_explanation('Unicycler now uses Racon to polish the miniasm assembly. It does '
                        'multiple rounds of polishing to get the best consensus, and circular '
                        'unitigs are rotated between rounds such that all parts (including the '
                        'ends) are polished well.')

    polish_dir = os.path.join(miniasm_dir, 'racon_polish')
    if not os.path.isdir(polish_dir):
        os.makedirs(polish_dir)

    # Save reads to file for polishing. We exclude reads that miniasm decided were chimeras.
    chimeric_read_names = set()
    chimeric_read_list_filename = os.path.join(miniasm_dir, 'chimeric_reads.txt')
    if os.path.isfile(chimeric_read_list_filename):
        with open(chimeric_read_list_filename) as chimeric_read_list_file:
            for line in chimeric_read_list_file:
                chimeric_read_names.add(line.strip())
    assembly_read_names = [x for x in assembly_read_names if x not in chimeric_read_names]
    polish_reads = os.path.join(polish_dir, 'polishing_reads.fastq')
    save_assembly_reads_to_file(polish_reads, assembly_read_names, read_dict, graph,
                                settings.RACON_CONTIG_DUPLICATION_COUNT)
    log.log('')

    col_widths = [6, 12, 14]
    racon_table_header = ['Polish round', 'Assembly size', 'Mapping quality']
    print_table([racon_table_header], fixed_col_widths=col_widths, left_align_header=False,
                alignments='RRR', indent=0)

    best_mapping_quality = 0
    best_fasta = None
    times_quality_failed_to_beat_best = 0
    for polish_round_count in range(settings.RACON_POLISH_LOOP_COUNT):
        round_num_str = '%02d' % (polish_round_count + 1)
        previous_round_num_str = '%02d' % polish_round_count
        pre_polish_fasta = os.path.join(polish_dir, round_num_str + '_1_pre-polish.fasta')
        mappings_filename = os.path.join(polish_dir, round_num_str + '_2_alignments.paf')
        racon_log = os.path.join(polish_dir, round_num_str + '_3_racon.log')
        post_polish_fasta = os.path.join(polish_dir, round_num_str + '_4_post-polish.fasta')
        if polish_round_count == 0:
            previous_fasta = pre_polish_fasta
        else:
            previous_fasta = os.path.join(polish_dir,
                                          previous_round_num_str + '_4_post-polish.fasta')
        if polish_round_count > 0:
            unitig_graph.rotate_circular_sequences()
        unitig_graph.save_to_fasta(pre_polish_fasta)

        # Create the alignments for Racon using miniasm.
        minimap_alignments_str = minimap_align_reads(pre_polish_fasta, polish_reads, 1, 0)
        alignments_by_read = load_minimap_alignments(minimap_alignments_str, filter_overlaps=True,
                                                     allowed_overlap=10, filter_by_minimisers=True)
        mapping_quality = 0
        read_alignments = load_minimap_alignments_basic(minimap_alignments_str)
        read_names = remove_dupes_preserve_order([x.read_name for x in read_alignments])
        with open(mappings_filename, 'wt') as mappings:
            for read_name in read_names:
                if read_name in alignments_by_read:
                    for a in alignments_by_read[read_name]:
                        mappings.write(a.paf_line)
                        mappings.write('\n')
                        mapping_quality += a.matching_bases

        racon_table_row = ['begin' if polish_round_count == 0 else str(polish_round_count),
                           int_to_str(unitig_graph.get_total_segment_length()),
                           int_to_str(mapping_quality)]
        print_table([racon_table_row], fixed_col_widths=col_widths, left_align_header=False,
                    alignments='LRR', indent=0, header_format='normal')


        # If we've failed to improve on our best quality for a few rounds, then we're done!
        if mapping_quality <= best_mapping_quality:
            times_quality_failed_to_beat_best += 1
        else:
            best_mapping_quality = mapping_quality
            best_fasta = previous_fasta
            times_quality_failed_to_beat_best = 0
        if times_quality_failed_to_beat_best > 2:
            break

        # Run Racon. It crashes sometimes, so repeat until its return code is 0.
        command = [racon_path, '--verbose', '9', '--threads', str(threads), '--bq', '-1',
                   polish_reads, mappings_filename, pre_polish_fasta, post_polish_fasta]
        return_code = 1
        for t in range(100):  # Only try a fixed number of times, to prevent an infinite loop.
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = process.communicate()
            with open(racon_log, 'wb') as log_file:
                log_file.write(out)
                log_file.write(err)
            return_code = process.returncode
            if return_code == 0 and os.path.isfile(post_polish_fasta):
                break
            if os.path.isfile(post_polish_fasta):
                os.remove(post_polish_fasta)
            if os.path.isfile(racon_log):
                os.remove(racon_log)

        # If even after all those tries Racon still didn't succeed, then we give up!
        if return_code != 0 or not os.path.isfile(post_polish_fasta):
            break

        unitig_graph.replace_with_polished_sequences(post_polish_fasta, scoring_scheme)

    log.log('')
    if best_fasta:
        log.log('Best polish: ' + best_fasta)
        unitig_graph.replace_with_polished_sequences(best_fasta, scoring_scheme)


def place_contigs(miniasm_dir, assembly_graph, unitig_graph, threads):
    log.log('', verbosity=2)
    log.log_explanation('Unicycler now places the single copy contigs back into the unitig '
                        'graph so it can extract bridging sequences between these contigs.',
                        verbosity=2)

    contig_search_dir = os.path.join(miniasm_dir, 'contig_search')
    if not os.path.isdir(contig_search_dir):
        os.makedirs(contig_search_dir)

    # Save the contigs to a FASTA file.
    contigs_fasta = os.path.join(contig_search_dir, 'contigs.fasta')
    contig_names = []
    longest_contig_seq_len = 0
    with open(contigs_fasta, 'wt') as fasta:
        for seg in sorted(assembly_graph.segments.values(), key=lambda x: x.number):
            if segment_suitable_for_miniasm_assembly(assembly_graph, seg):
                contig_name = 'CONTIG_' + str(seg.number)
                contig_names.append(contig_name)
                fasta.write('>' + contig_name + '\n')
                seq = seg.forward_sequence
                fasta.write(seq)
                longest_contig_seq_len = max(longest_contig_seq_len, len(seq))
                fasta.write('\n')

    # Save the unitigs to a FASTA file. If the segment is circular, then we save it extra sequence
    # on the end, so we can find contigs which loop. We also save these extended sequences in a
    # dictionary for later.
    unitigs_fasta = os.path.join(contig_search_dir, 'unitigs.fasta')
    extended_unitig_sequences = {}
    with open(unitigs_fasta, 'wt') as fasta:
        for seg in sorted(unitig_graph.segments.values(), key=lambda x: x.get_length(),
                          reverse=True):
            seg_seq = seg.forward_sequence
            if unitig_graph.segment_is_circular(seg.full_name):
                if len(seg_seq) <= longest_contig_seq_len:
                    seg_seq += seg_seq
                else:
                    seg_seq += seg_seq[:longest_contig_seq_len]
            fasta.write('>' + seg.full_name + '\n')
            fasta.write(seg_seq)
            fasta.write('\n')
            extended_unitig_sequences[seg.full_name] = seg_seq

    alignments = load_minimap_alignments(minimap_align_reads(unitigs_fasta, contigs_fasta,
                                                             threads, 0, 'find contigs'))
    for contig_name in contig_names:
        contig_alignments = alignments[contig_name]

        # Clean up the alignments by removing redundant/conflicting hits and combining separated
        # pieces.
        filtered_alignments = []
        for a in contig_alignments:
            unitig_name = a.ref_name
            unitig_length = unitig_graph.segments[unitig_name].get_length()
            if unitig_graph.segment_is_circular(unitig_name) and a.ref_end > unitig_length:
                a.ref_start -= unitig_length
                a.ref_end -= unitig_length
            filtered_alignments.append(a)
        contig_alignments = remove_conflicting_alignments(filtered_alignments, 10)
        contig_alignments = combine_close_hits(contig_alignments, settings.FOUND_CONTIG_MIN_RATIO,
                                               settings.FOUND_CONTIG_MAX_RATIO)
        if not contig_alignments:
            continue
        a = sorted(contig_alignments, key=lambda x: x.fraction_read_aligned(), reverse=True)[0]
        if a.fraction_read_aligned() < settings.FOUND_CONTIG_MIN_FRACTION:
            continue
        start_gap = a.read_start
        end_gap = a.read_length - a.read_end
        estimated_unitig_start_pos = a.ref_start - start_gap
        estimated_unitig_end_pos = a.ref_end + end_gap
        if start_gap > settings.FOUND_CONTIG_MAX_END_GAP or \
                end_gap > settings.FOUND_CONTIG_MAX_END_GAP:
            continue
        unitig_name = a.ref_name
        if estimated_unitig_start_pos < 0:
            if unitig_graph.segment_is_circular(unitig_name):
                a.ref_start += unitig_length
                a.ref_end += unitig_length
                estimated_unitig_start_pos += unitig_length
                estimated_unitig_end_pos += unitig_length
            else:
                estimated_unitig_start_pos = 0


        print(contig_name + ': ' + str(a))  # TEMP

        # Now we have one and only one alignment for this contig, but it's a minimap alignment that
        # is likely trimming a bit off the start/end. We now extend the alignment to cover the
        # entire contig.
        extended_unitig_seq = extended_unitig_sequences[a.ref_name]
        gap = max(settings.FOUND_CONTIG_END_REFINEMENT_SIZE, 2 * max(start_gap, end_gap))
        contig_number = int(a.read_name[7:])
        contig_start_seq = assembly_graph.segments[contig_number].forward_sequence[:gap]
        contig_end_seq = assembly_graph.segments[contig_number].forward_sequence[-gap:]
        if a.read_strand == '-':
            start_gap, end_gap = end_gap, start_gap
            contig_start_seq, contig_end_seq = \
                reverse_complement(contig_end_seq), reverse_complement(contig_start_seq)

        start_region_start = max(0, estimated_unitig_start_pos - gap)
        start_region_end = min(len(extended_unitig_seq), estimated_unitig_start_pos + gap)
        unitig_start_region = extended_unitig_seq[start_region_start:start_region_end]

        end_region_start = max(0, estimated_unitig_end_pos - gap)
        end_region_end = min(len(extended_unitig_seq), estimated_unitig_end_pos + gap)
        unitig_end_region = extended_unitig_seq[end_region_start:end_region_end]

        print('start_gap:', start_gap)  # TEMP
        print('end_gap:', end_gap)  # TEMP
        print('contig_start_seq:', contig_start_seq)  # TEMP
        print('contig_end_seq:', contig_end_seq)  # TEMP
        print('estimated_unitig_start_pos:', estimated_unitig_start_pos)  # TEMP
        print('estimated_unitig_end_pos:', estimated_unitig_end_pos)  # TEMP
        print('unitig_start_region:', unitig_start_region)  # TEMP
        print('unitig_end_region:', unitig_end_region)  # TEMP

        print('')  # TEMP





