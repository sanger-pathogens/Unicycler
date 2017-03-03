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

from . import log


def build_miniasm_bridges(graph, out_dir, keep, threads, read_dict, long_read_filename):

    log.log_section_header('Assemble contigs and long reads with miniasm and Racon')
    log.log_explanation('Unicycler uses miniasm and Racon to construct a string graph assembly '
                        'using both the Illumina contigs and the long reads. If this produces a '
                        'reliable assembly, Unicycler can extract bridges between contigs and use '
                        'them to simplify the assembly graph.')
    log.log_explanation('This method requires decent coverage of long reads and may therefore not '
                        'be fruitful if long reads are sparse. However, this method does not rely '
                        'on the short read assembly graph having good connectivity and is '
                        'therefore able to bridge an assembly graph even when it contains many '
                        'dead ends.')

    # ALIGN LONG READS TO THE SINGLE COPY CONTIGS USING MINIMAP.

    # EXTRACT READS USEFUL FOR LONG READ ASSEMBLY.
    # * Take all single copy contigs over a certain length and get reads which overlap two or more.
    #   * Some logic is currently in get_overlapping_reads.py.
    #   * While I'm at it, I should throw out reads which look like chimeras based on incompatible
    #     mapping.
    # * Create a file of "long reads" which contains:
    #   * real long reads as found (and possibly split) by the above step
    #   * single copy contigs in FASTQ form (with a high quality, 'I' or something)

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