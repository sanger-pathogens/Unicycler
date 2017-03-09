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
from .minimap_alignment import align_long_reads_to_assembly_graph, load_minimap_alignments_basic, \
    build_start_end_overlap_sets
from .cpp_wrappers import minimap_align_reads, miniasm_assembly
from .string_graph import StringGraph
from . import log
from . import settings


class MiniasmFailure(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


def build_miniasm_bridges(graph, out_dir, keep, threads, read_dict, long_read_filename):
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
                        'may not be fruitful if long reads are sparse. However, this method does '
                        'not rely on the short read assembly graph having good connectivity and is '
                        'able to bridge an assembly graph even when it contains many dead ends.',
                        extra_empty_lines_after=0)
    log.log_explanation('Unicycler uses two types of "reads" as assembly input: '
                        'sufficiently long single-copy short read contigs and actual long reads '
                        'which overlap two or more of these contigs. It then assembles them with '
                        'a modified version of miniasm which gives precedence to the contigs over '
                        'the real long reads.')

    miniasm_dir = os.path.join(out_dir, 'miniasm_assembly')
    if not os.path.exists(miniasm_dir):
        os.makedirs(miniasm_dir)

    # Align all the long reads to the graph and get the ones which overlap single-copy contigs
    # (and will therefore be useful for assembly).
    minimap_alignments = align_long_reads_to_assembly_graph(graph, long_read_filename,
                                                            miniasm_dir, threads)
    start_overlap_reads, end_overlap_reads = build_start_end_overlap_sets(minimap_alignments)
    assembly_read_names = get_miniasm_assembly_reads(minimap_alignments, graph)

    assembly_reads_filename = os.path.join(miniasm_dir, '01_assembly_reads.fastq')
    mappings_filename = os.path.join(miniasm_dir, '02_mappings.paf')

    # Look at the alignments to see where we need to split the reads to ensure no contigs are
    # contained in reads.
    read_split_points = {}
    for read_name, alignments in minimap_alignments.items():
        for a in alignments:
            seg = graph.segments[int(a.ref_name)]
            if segment_suitable_for_miniasm_assembly(graph, seg) and a.ref_contained_in_read():
                split_start = a.read_start + settings.BROKEN_ASSEMBLY_READ_END_GAP
                split_end = a.read_end - settings.BROKEN_ASSEMBLY_READ_END_GAP
                if read_name not in read_split_points:
                    read_split_points[read_name] = []
                read_split_points[read_name].append((split_start, split_end))

    while True:
        # Prepare an assembly FASTQ which contains both contigs and real long reads. The real long
        # reads may be split to ensure that no contigs are contained within read (or else miniasm
        # will filter them out of its assembly).
        mean_read_quals = save_assembly_reads_to_file(assembly_reads_filename, assembly_read_names,
                                                      read_dict, graph, read_split_points)

        # Do an all-vs-all alignment of the assembly FASTQ, for miniasm input.
        minimap_alignments_str = minimap_align_reads(assembly_reads_filename,
                                                     assembly_reads_filename, threads, 0, True)

        # TO DO: FILTER OUT ALIGNMENTS WITH WEIRD OVERHANGS ON THE CONTIGS
        # TO DO: FILTER OUT ALIGNMENTS WITH WEIRD OVERHANGS ON THE CONTIGS
        # TO DO: FILTER OUT ALIGNMENTS WITH WEIRD OVERHANGS ON THE CONTIGS
        # TO DO: FILTER OUT ALIGNMENTS WITH WEIRD OVERHANGS ON THE CONTIGS
        # TO DO: FILTER OUT ALIGNMENTS WITH WEIRD OVERHANGS ON THE CONTIGS
        # TO DO: FILTER OUT ALIGNMENTS WITH WEIRD OVERHANGS ON THE CONTIGS

        # Before we continue, we need to double check that no contigs are contained within reads in
        # our new alignment. If they are, add new split points and repeat the above process.
        more_read_split_points = contigs_are_contained_in_long_reads(minimap_alignments_str)
        if more_read_split_points:
            read_split_points = add_more_read_split_points(read_split_points,
                                                           more_read_split_points)

        # If no contigs are contained, then we are good to continue!
        else:
            with open(mappings_filename, 'wt') as mappings:
                mappings.write(minimap_alignments_str)
            break

    # Now actually do the miniasm assembly, which will create a GFA file of the string graph.
    # TO DO: intelligently set the min_ovlp setting (currently 1) based on the depth? The miniasm
    #        default is 3, so perhaps use 3 if the depth is high enough, 2 if it's lower and 1 if
    #        it's very low. I'm not yet sure what the risks are (if any) with using a min_ovlp of
    #        1 when the depth is high.
    log.log('Assembling reads with miniasm')
    miniasm_assembly(assembly_reads_filename, mappings_filename, miniasm_dir)
    before_transitive_reduction_filename = os.path.join(miniasm_dir, '03_raw_string_graph.gfa')
    string_graph_filename = os.path.join(miniasm_dir, '10_final_string_graph.gfa')
    if not (os.path.isfile(string_graph_filename) and
            os.path.isfile(before_transitive_reduction_filename)):
        raise MiniasmFailure('miniasm failed to generate a string graph')
    string_graph = StringGraph(string_graph_filename, mean_read_quals)
    before_transitive_reduction = StringGraph(before_transitive_reduction_filename, mean_read_quals)

    string_graph.remove_non_bridging_paths()
    string_graph.save_to_gfa(os.path.join(miniasm_dir, '12_non_bridging_paths_removed.gfa'))
    string_graph.remove_overlaps(before_transitive_reduction)
    string_graph.save_to_gfa(os.path.join(miniasm_dir, '13_overlaps_removed.gfa'))
    string_graph.merge_reads()
    string_graph.save_to_gfa(os.path.join(miniasm_dir, '14_reads_merged.gfa'))







    # EXTRACT ALL CONTIG-CONTIG BRIDGE SEQUENCES.
    # * Any two single copy contigs connected by an unbranching path that contains no other contigs.

    # POLISH EACH BRIDGE SEQUENCE.
    # * For this we use the set of long reads which overlap the two single copy contigs on the
    #   correct side. It is not necessary for reads to overlap both contigs, as this will give us
    #   better coverage in the intervening repeat region.
    # * Use only the long read sequences, not the Illumina contigs. Since the Illumina contigs may
    #   not have been used all the way to their ends (slightly trimmed), this means a bit of contig
    #   sequence may be replaced by long read consensus.

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


def contigs_are_contained_in_long_reads(minimap_alignments_str):
    """
    This function takes the minimap output of an all-vs-all alignment for miniasm. If it finds that
    any contig is contained within a read, it returns them in a list.
    """
    more_read_split_points = {}
    minimap_alignments = load_minimap_alignments_basic(minimap_alignments_str)

    for a in minimap_alignments:

        # Since this alignment is all-vs-all, either the 'read' or the 'ref' could be a contig.
        ref_is_contig = a.ref_name.startswith('CONTIG_')
        read_is_contig = a.read_name.startswith('CONTIG_')

        # This covers where the 'read' is an actual read and the 'ref' is a contig.
        if ref_is_contig and not read_is_contig and a.ref_contained_in_read():
            subread_start_point = int(a.read_name.split(':')[-1].split('-')[0])
            split_start = a.read_start + subread_start_point + settings.BROKEN_ASSEMBLY_READ_END_GAP
            split_end = a.read_end + subread_start_point - settings.BROKEN_ASSEMBLY_READ_END_GAP
            if a.read_name not in more_read_split_points:
                more_read_split_points[a.read_name] = []
            more_read_split_points[a.read_name].append((split_start, split_end))

        # This covers where the 'read' is a contig and the 'ref' is an actual read.
        if read_is_contig and not ref_is_contig and a.read_contained_in_ref():
            subread_start_point = int(a.ref_name.split(':')[-1].split('-')[0])
            split_start = a.ref_start + subread_start_point + settings.BROKEN_ASSEMBLY_READ_END_GAP
            split_end = a.ref_end + subread_start_point - settings.BROKEN_ASSEMBLY_READ_END_GAP
            if a.ref_name not in more_read_split_points:
                more_read_split_points[a.ref_name] = []
            more_read_split_points[a.ref_name].append((split_start, split_end))

    return more_read_split_points


def add_more_read_split_points(read_split_points, new_split_points):
    for read_name, split_points in new_split_points.items():
        read_name_no_range = read_name.rsplit(':', 1)[0]
        if read_name_no_range not in read_split_points:
            read_split_points[read_name_no_range] = []
        read_split_points[read_name_no_range] = sorted(read_split_points[read_name_no_range] +
                                                       split_points)

    # Now we process the split points to deal with overlaps and containments.
    for read_name, split_points in read_split_points.items():
        points = []
        for split_point in split_points:
            points.append((split_point[0], 1))
            points.append((split_point[1], -1))
        points = sorted(points)

        new_points = []
        new_start = None
        for p, change in points:
            if change == 1:
                new_start = p
            else:  # change == -1:
                if new_start is not None:
                    new_points.append((new_start, p))
                new_start = None
        read_split_points[read_name] = new_points

    return read_split_points


def get_miniasm_assembly_reads(minimap_alignments, graph):
    """
    Returns a list of read names which overlap at least two different single copy graph segments.
    """
    miniasm_assembly_reads = []
    for read_name, alignments in minimap_alignments.items():
        overlap_count = 0
        for a in alignments:
            if a.overlaps_reference():
                seg = graph.segments[int(a.ref_name)]
                if segment_suitable_for_miniasm_assembly(graph, seg):
                    overlap_count += 1

        # TO DO: I'm not sure if this value should be 2 (only taking reads which span all the way
        # from one contig to the next) or 1 (also taking reads which overlap one contig but do not
        # reach the next). I should test each option and analyse.
        if overlap_count >= 1:
            miniasm_assembly_reads.append(read_name)
    return sorted(miniasm_assembly_reads)


def save_assembly_reads_to_file(read_filename, read_names, read_dict, graph, read_split_points):
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

        # Now save the actual long reads (though they may be split to prevent Illumina contigs
        # from being contained).
        for read_name in read_names:
            read = read_dict[read_name]
            read_length = read.get_length()
            if read_name not in read_split_points:
                ranges = [(0, read_length)]
            else:
                range_starts, range_ends = [], []
                for split_point in read_split_points[read_name]:
                    range_starts.append(split_point[0])
                    range_ends.append(split_point[1])
                range_starts = [0] + range_starts
                range_ends.append(read_length)
                ranges = list(zip(range_starts, range_ends))

            for i, out_range in enumerate(ranges):
                s, e = out_range[0], out_range[1]
                read_part_name = read_name + ':' + str(s) + '-' + str(e)
                seq = read.sequence[s:e]
                quals = read.qualities[s:e]
                fastq.write('@')
                fastq.write(read_part_name)
                fastq.write('\n')
                fastq.write(seq)
                fastq.write('\n+\n')
                fastq.write(quals)
                fastq.write('\n')
                mean_read_quals[read_part_name] = statistics.mean(ord(x)-33 for x in quals)
        log.log('  ' + str(len(read_names)) + ' overlapping long reads (out of ' +
                str(len(read_dict)) + ' total long reads)')
    log.log('')
    return mean_read_quals


# def get_assembly_output_ranges(read_alignments, read_length, graph):
#     """
#     This function outputs the part(s) of the read which should be output as reads for assembly.
#     The whole read isn't used because we don't want the read to contain a graph segment, because
#     that would result in miniasm throwing out the graph segment because it is contained.
#     """
#     range_starts, range_ends = [], []
#     for a in read_alignments:
#         seg = graph.segments[int(a.ref_name)]
#         if segment_suitable_for_miniasm_assembly(graph, seg) and a.ref_contained_in_read():
#             range_starts.append(a.read_start + settings.BROKEN_ASSEMBLY_READ_END_GAP)
#             range_ends.append(a.read_end - settings.BROKEN_ASSEMBLY_READ_END_GAP)
#     range_starts = [0] + range_starts
#     range_ends.append(read_length)
#     return list(zip(range_starts, range_ends))


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
