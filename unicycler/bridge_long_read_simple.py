"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module takes care of simple long read bridging. This is where parts of the graph have a very
straightforward repeat and minimap alignments can conclusively support one over the alternatives.
Simple long read bridging is done before other, more complex bridging processes.

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
from collections import defaultdict
import itertools
from .cpp_wrappers import minimap_align_reads
from .minimap_alignment import load_minimap_alignments
from .read_ref import load_references
from . import log
from . import settings


def apply_simple_long_read_bridges(graph, out_dir, keep, threads, read_dict, read_names,
                                   long_read_filename):
    """
    LOOK FOR VERY SIMPLE PARTS IN THE GRAPH WHICH CAN BE UNAMBIGUOUSLY RESOLVED WITH THE LONG
    READS.
    * For example, a two-way merge and split.
    * Places in the graph where there are X ways in, X ways out, and each way in to way out has
      one and only one path.
    * First, find these locations.
    * Second, determine all possible paths. For a 2-2 junction, this will be 4 paths. For a 3-3
      junction, this will be 9 paths, etc.
    * Evaluate each path using the minimap alignments and sort the paths from best to worst.
    * If the top X paths (where X is the in and out degree) are complete and not conflicting
      (they use each in path once and each out path once), then we can apply the bridges!
    """
    log.log_section_header('Simple repeat bridging with minimap alignments')

    bridging_dir = os.path.join(out_dir, 'simple_bridging')
    if not os.path.exists(bridging_dir):
        os.makedirs(bridging_dir)

    # Do a minimap alignment of all reads to all segments.
    segments_fasta = os.path.join(bridging_dir, 'all_segments.fasta')
    graph.save_to_fasta(segments_fasta)
    references = load_references(segments_fasta, section_header=None)
    reference_dict = {x.name: x for x in references}
    minimap_alignments_str = minimap_align_reads(segments_fasta, long_read_filename, threads, 0)
    minimap_alignments = load_minimap_alignments(minimap_alignments_str, read_dict, reference_dict,
                                                 filter_overlaps=True,
                                                 allowed_overlap=settings.ALLOWED_MINIMAP_OVERLAP,
                                                 filter_by_minimisers=True,
                                                 minimiser_ratio=settings.MAX_TO_MIN_MINIMISER_RATIO)

    # Build indices of start and end contig overlaps so we can quickly determine which reads overlap
    # which end of a particular contig.
    # These dictionaries have a key of a signed segment number and a value of a set of reads names.
    start_overlap_reads = defaultdict(set)
    end_overlap_reads = defaultdict(set)
    min_overlap_amount = 100
    for read_name, alignments in minimap_alignments.items():
        for a in alignments:
            seg_num = int(a.ref_name)
            # print()  # TEMP
            # print(read_name + ' (' + str(a.read_length) + ' bp)')  # TEMP
            # print(a)  # TEMP
            if a.read_strand == '+':
                seg_start = a.ref_start
                seg_end = a.ref_end
            else:  # a.read_strand == '-'
                seg_num *= -1
                seg_start = a.ref_length - a.ref_end
                seg_end = a.ref_length - a.ref_start
            # print('seg_start =', seg_start)  # TEMP
            # print('seg_end =', seg_end)  # TEMP
            adjusted_seg_start = seg_start - a.read_start
            adjusted_seg_end = seg_end + a.read_end_gap
            # print('adjusted_seg_start =', adjusted_seg_start)  # TEMP
            # print('adjusted_seg_end =', adjusted_seg_end)  # TEMP
            if adjusted_seg_start < -min_overlap_amount:
                start_overlap_reads[seg_num].add(read_name)
            if adjusted_seg_end > a.ref_length + min_overlap_amount:
                end_overlap_reads[seg_num].add(read_name)

    # Simple bridge two-way junctions
    junctions = graph.find_simple_two_way_junctions(settings.MIN_SEGMENT_LENGTH_FOR_SIMPLE_BRIDGING)
    for junction in junctions:
        inputs = graph.reverse_links[junction]
        outputs = graph.forward_links[junction]

        # There are two possible ways to resolve a simple two-way junction. Determine the paths for
        # each option.
        option_1_paths = [[inputs[0], junction, outputs[0]],
                          [inputs[1], junction, outputs[1]],
                          [-outputs[0], -junction, -inputs[0]],
                          [-outputs[1], -junction, -inputs[1]]]
        option_2_paths = [[inputs[0], junction, outputs[1]],
                          [inputs[1], junction, outputs[0]],
                          [-outputs[0], -junction, -inputs[1]],
                          [-outputs[1], -junction, -inputs[0]]]

        # # Make a FASTA file of the four possible paths (only need one strand).
        # options_fasta_filename = os.path.join(bridging_dir,
        #                                       'junction_' + str(junction) + '_options.fasta')
        # with open(options_fasta_filename, 'wt') as options_fasta:
        #     options_fasta.write('>option_1_path_1\n')
        #     options_fasta.write(graph.get_path_sequence(option_1_paths[0]) + '\n')
        #     options_fasta.write('>option_1_path_2\n')
        #     options_fasta.write(graph.get_path_sequence(option_1_paths[1]) + '\n')
        #     options_fasta.write('>option_2_path_1\n')
        #     options_fasta.write(graph.get_path_sequence(option_2_paths[0]) + '\n')
        #     options_fasta.write('>option_2_path_2\n')
        #     options_fasta.write(graph.get_path_sequence(option_2_paths[1]) + '\n')

        # Gather up all reads which could possibly span this junction.
        relevant_reads = list(end_overlap_reads[inputs[0]] | end_overlap_reads[inputs[1]] |
                              end_overlap_reads[-outputs[0]] | end_overlap_reads[-outputs[1]] |
                              start_overlap_reads[outputs[0]] | start_overlap_reads[outputs[1]] |
                              start_overlap_reads[-inputs[0]] | start_overlap_reads[-inputs[1]])


        print()  # TEMP
        print()  # TEMP
        print()  # TEMP
        print()  # TEMP
        print()  # TEMP
        print(junction)  # TEMP
        print()  # TEMP


        # Each read now casts a vote in one of four ways:
        #   1) Option 1: the read supports the connection of option 1
        #   2) Option 2: the read supports the connection of option 2
        #   3) Neither option: the read supports some other strange connection
        #   4) No vote: the read doesn't support any connection
        option_1_votes = 0
        option_2_votes = 0
        neither_option_votes = 0

        expected_next_seg = [[inputs[0], outputs[0], outputs[1]],
                             [inputs[1], outputs[1], outputs[0]],
                             [-outputs[0], -inputs[0], -inputs[1]],
                             [-outputs[1], -inputs[1], -inputs[0]]]
        for r in relevant_reads:
            # print(r, end='  ')  # TEMP

            # Turn the alignments into a list of segment numbers, excluding the junction segment.
            alignments = [int(x.ref_name) * (-1 if x.read_strand == '-' else 1)
                          for x in minimap_alignments[r] if x.ref_name != str(junction)]

            # Remove adjacent duplicates in the alignments:
            alignments = [k for k, g in itertools.groupby(alignments)]
            # print(alignments)  # TEMP

            for p in expected_next_seg:
                start, option_1_end, option_2_end = p[0], p[1], p[2]
                try:
                    after_start = alignments[alignments.index(start) + 1]
                    if after_start == option_1_end:
                        option_1_votes += 1
                        # print('VOTE: 1  ', start, option_1_end)  # TEMP
                    elif after_start == option_2_end:
                        option_2_votes += 1
                        # print('VOTE: 2  ', start, option_2_end)  # TEMP
                    else:
                        neither_option_votes += 1
                        print(r)  # TEMP
                        for a in minimap_alignments[r]:  # TEMP
                            print(a)  # TEMP
                        print('VOTE: NEITHER')  # TEMP
                except (ValueError, IndexError):
                    pass

        print()  # TEMP
        print('option_1_votes:', option_1_votes)  # TEMP
        print('option_2_votes:', option_2_votes)  # TEMP
        print('neither_option_votes:', neither_option_votes)  # TEMP
        print()  # TEMP











    if keep < 3:
        shutil.rmtree(bridging_dir)


