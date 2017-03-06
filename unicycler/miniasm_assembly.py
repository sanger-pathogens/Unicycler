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
from .minimap_alignment import align_long_reads_to_assembly_graph
from .cpp_wrappers import minimap_align_reads, miniasm_assembly
from . import log
from . import settings


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

    minimap_alignments = align_long_reads_to_assembly_graph(graph, long_read_filename, miniasm_dir,
                                                            threads, read_dict)
    assembly_read_names = get_miniasm_assembly_reads(minimap_alignments, graph)

    # Save appropriate single copy contigs and informative reads to a FASTQ file.
    long_read_filename = os.path.join(miniasm_dir, '01_assembly_reads.fastq')
    save_assembly_reads_to_file(minimap_alignments, long_read_filename, assembly_read_names,
                                read_dict, graph)

    # Do an all vs all alignment and save the results in a PAF file.
    log.log('Finding read-read overlaps using minimap')
    minimap_alignments_str = minimap_align_reads(long_read_filename, long_read_filename, threads,
                                                 0, True)
    mappings_filename = os.path.join(miniasm_dir, '02_mappings.paf')
    log.log('Saving ' + mappings_filename)
    with open(mappings_filename, 'wt') as mappings:
        mappings.write(minimap_alignments_str)

    # Now actually do the miniasm assembly, which will create a GFA file of the string graph.
    log.log('Assembling reads with miniasm')
    miniasm_assembly(long_read_filename, mappings_filename, miniasm_dir)

    # ASSEMBLE LONG READS USING MINIASM.
    # * The min_ovlp setting should either be 1 or dynamically determined based on depth.
    # * Final target: the string graph.

    # REMOVE OVERLAPS FROM MINIASM STRING GRAPH.
    # * Selectively remove overlaps from lower quality sequences first.
    # * Keep as much contig sequence as possible.
    # * Process idea:
    #   * Find the lowest quality read, based on average qscore and remove as much as possible
    #     (could be all of the read if its neighbours overlap).
    #   * Repeat until there are no more overlaps.
    # * Merge the assembly together (keeping single copy contigs separate).

    # EXTRACT ALL CONTIG-CONTIG BRIDGE SEQUENCES.
    # * Any two single copy contigs connected by an unbranching path that contains no other contigs.

    # POLISH EACH BRIDGE SEQUENCE.
    # * For this we use the set of long reads which overlap the two single copy contigs on the
    #   correct side. It is not necessary for reads to overlap both contigs, as this will give us
    #   better coverage in the intervening repeat region.
    # * Include the single copy contigs as 'reads'.
    #   * Specifically, use the slightly trimmed single copy contigs in the miniasm string graph
    #     (in case the entire contig has a bogus end).
    #   * The high qscores here should ensure that no changes are made in these regions. Include
    #     extra copies if necessary.

    # MAKE EACH BRIDGE SEGMENT.
    # * Goal is to turn one contiguous sequence spanning both single copy contigs into:
    #   CONTIG -> BRIDGE -> CONTIG

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
        if overlap_count >= 2:
            miniasm_assembly_reads.append(read_name)
    return sorted(miniasm_assembly_reads)


def save_assembly_reads_to_file(minimap_alignments, read_filename, read_names, read_dict, graph):
    qual = chr(settings.CONTIG_READ_QSCORE + 33)
    log.log('Saving to ' + read_filename + ':')

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
            ranges = get_assembly_output_ranges(minimap_alignments[read_name], read.get_length(),
                                                graph)
            for i, out_range in enumerate(ranges):
                s, e = out_range[0], out_range[1]
                fastq.write('@')
                fastq.write(read_name + '_' + str(i))
                fastq.write('\n')
                fastq.write(read.sequence[s:e])
                fastq.write('\n+\n')
                fastq.write(read.qualities[s:e])
                fastq.write('\n')
        log.log('  ' + str(len(read_names)) + ' overlapping long reads (out of ' +
                str(len(read_dict)) + ' total long reads)')
    log.log('')


def get_assembly_output_ranges(read_alignments, read_length, graph):
    """
    This function outputs the part(s) of the read which should be output as reads for assembly.
    The whole read isn't used because we don't want the read to contain a graph segment, because
    that would result in miniasm throwing out the graph segment because it is contained.
    """
    range_starts, range_ends = [], []
    for a in read_alignments:
        seg = graph.segments[int(a.ref_name)]
        if segment_suitable_for_miniasm_assembly(graph, seg) and a.ref_contained_in_read():
            range_starts.append(a.read_start + settings.BROKEN_ASSEMBLY_READ_END_GAP)
            range_ends.append(a.read_end - settings.BROKEN_ASSEMBLY_READ_END_GAP)
    range_starts = [0] + range_starts
    range_ends.append(read_length)
    return list(zip(range_starts, range_ends))


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
