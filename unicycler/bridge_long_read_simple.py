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
from .cpp_wrappers import minimap_align_reads, fully_global_alignment
from .minimap_alignment import load_minimap_alignments
from .read_ref import load_references
from .assembly_graph_segment import Segment
from .assembly_graph_copy_depth import determine_copy_depth
from .misc import print_table, get_right_arrow
from .alignment import AlignmentScoringScheme
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
    log.log('Aligning long reads to graph using minimap', 2)
    graph.save_to_fasta(segments_fasta, verbosity=2)
    references = load_references(segments_fasta, section_header=None, show_progress=False)
    reference_dict = {x.name: x for x in references}
    minimap_alignments_str = minimap_align_reads(segments_fasta, long_read_filename, threads, 0)
    minimap_alignments = load_minimap_alignments(minimap_alignments_str, read_dict, reference_dict,
                                                 filter_overlaps=True,
                                                 allowed_overlap=settings.ALLOWED_MINIMAP_OVERLAP,
                                                 filter_by_minimisers=True,
                                                 minimiser_ratio=settings.MAX_TO_MIN_MINIMISER_RATIO)
    log.log('Number of minimap alignments: ' + str(len(minimap_alignments)), 2)
    log.log('', 2)

    # Build indices of start and end contig overlaps so we can quickly determine which reads
    # overlap which end of a particular contig.
    # These dictionaries have a key of a signed segment number and a value of a set of reads names.
    start_overlap_reads = defaultdict(set)
    end_overlap_reads = defaultdict(set)
    min_overlap_amount = 100
    for read_name, alignments in minimap_alignments.items():
        for a in alignments:
            seg_num = int(a.ref_name)
            if a.read_strand == '+':
                seg_start = a.ref_start
                seg_end = a.ref_end
            else:  # a.read_strand == '-'
                seg_num *= -1
                seg_start = a.ref_length - a.ref_end
                seg_end = a.ref_length - a.ref_start
            adjusted_seg_start = seg_start - a.read_start
            adjusted_seg_end = seg_end + a.read_end_gap
            if adjusted_seg_start < -min_overlap_amount:
                start_overlap_reads[seg_num].add(read_name)
            if adjusted_seg_end > a.ref_length + min_overlap_amount:
                end_overlap_reads[seg_num].add(read_name)

    simple_bridge_two_way_junctions(graph, start_overlap_reads, end_overlap_reads,
                                    minimap_alignments)
    simple_bridge_loops(graph, start_overlap_reads, end_overlap_reads, minimap_alignments,
                        read_dict)
    log.log('', 2)

    graph.merge_all_possible(None, 2)
    determine_copy_depth(graph)

    if keep < 3:
        shutil.rmtree(bridging_dir)
    return graph


def simple_bridge_two_way_junctions(graph, start_overlap_reads, end_overlap_reads,
                                    minimap_alignments):
    two_way_junctions_table = [['Junction', 'Option 1', 'Option 2', '1 votes', '2 votes',
                                'Alt votes', 'Result']]
    log.log('Simple two-way junction bridges:')

    junctions = graph.find_simple_two_way_junctions(settings.MIN_SEGMENT_LENGTH_FOR_SIMPLE_BRIDGING)
    for junction in junctions:
        table_row = [junction]

        inputs = graph.reverse_links[junction]
        outputs = graph.forward_links[junction]

        # There are two possible ways to resolve a simple two-way junction:
        #   1) inputs[0], junction, outputs[0]
        #      inputs[1], junction, outputs[1]
        #   2) inputs[0], junction, outputs[1]
        #      inputs[1], junction, outputs[0]]

        option_1_1_str = get_right_arrow().join(str(x) for x in [inputs[0], junction, outputs[0]])
        option_1_2_str = get_right_arrow().join(str(x) for x in [inputs[1], junction, outputs[1]])
        table_row.append(option_1_1_str + ' ' + option_1_2_str)
        option_2_1_str = get_right_arrow().join(str(x) for x in [inputs[0], junction, outputs[1]])
        option_2_2_str = get_right_arrow().join(str(x) for x in [inputs[1], junction, outputs[0]])
        table_row.append(option_2_1_str + ' ' + option_2_2_str)

        # Gather up all reads which could possibly span this junction.
        relevant_reads = list(end_overlap_reads[inputs[0]] | end_overlap_reads[inputs[1]] |
                              end_overlap_reads[-outputs[0]] | end_overlap_reads[-outputs[1]] |
                              start_overlap_reads[outputs[0]] | start_overlap_reads[outputs[1]] |
                              start_overlap_reads[-inputs[0]] | start_overlap_reads[-inputs[1]])

        # Each read now casts a vote in one of four ways:
        #   1) Option 1: the read supports the connection of option 1
        #   2) Option 2: the read supports the connection of option 2
        #   3) Neither option: the read supports some other strange connection
        #   4) No vote: the read doesn't support any connection
        option_1_votes = 0
        option_2_votes = 0
        neither_option_votes = 0

        # The lists in expected_next_seg each hold three segment numbers:
        #   1) The starting segment
        #   2) The segment that should follow if option 1 is correct.
        #   3) The segment that should follow if option 2 is correct.
        expected_next_seg = [[inputs[0], outputs[0], outputs[1]],
                             [inputs[1], outputs[1], outputs[0]],
                             [-outputs[0], -inputs[0], -inputs[1]],
                             [-outputs[1], -inputs[1], -inputs[0]]]

        for r in relevant_reads:

            # Turn the alignments into a list of segment numbers, excluding the junction segment.
            alignments = [int(x.ref_name) * (-1 if x.read_strand == '-' else 1)
                          for x in minimap_alignments[r] if x.ref_name != str(junction)]
            alignments = [k for k, g in itertools.groupby(alignments)]  # remove adjacent duplicates

            for p in expected_next_seg:
                start, option_1_end, option_2_end = p[0], p[1], p[2]
                try:
                    after_start = alignments[alignments.index(start) + 1]
                    if after_start == option_1_end:
                        option_1_votes += 1
                    elif after_start == option_2_end:
                        option_2_votes += 1
                    else:
                        neither_option_votes += 1
                except (ValueError, IndexError):
                    pass

        table_row += [str(option_1_votes), str(option_2_votes), str(neither_option_votes)]

        # If the neither option got the most votes, then something is screwy (e.g. perhaps this
        # isn't the simple two-way junction we thought), so we don't bridge anything.
        if neither_option_votes > option_1_votes and neither_option_votes > option_2_votes:
            table_row.append('too many alt')

        # If there aren't any votes at all, that's a shame! The long reads were probably too few
        # and/or too short to span the junction.
        elif option_1_votes == 0 and option_2_votes == 0:
            table_row.append('no support')

        # If the best option isn't sufficiently better than the second best option, we don't bridge
        # anything.
        elif min(option_1_votes, option_2_votes) / max(option_1_votes, option_2_votes) > 0.5:
            table_row.append('too close')

        else:
            # If we got here, then we're good to bridge!
            junction_seg = graph.segments[junction]
            bridge_depth = junction_seg.depth / 2
            bridge_seq = junction_seg.forward_sequence
            graph.remove_segments([junction])

            junction_copy_num_1 = graph.get_next_available_seg_number()
            junction_copy_num_2 = junction_copy_num_1 + 1
            junction_copy_1 = Segment(junction_copy_num_1, bridge_depth, bridge_seq, True)
            junction_copy_2 = Segment(junction_copy_num_2, bridge_depth, bridge_seq, True)
            graph.segments[junction_copy_num_1] = junction_copy_1
            graph.segments[junction_copy_num_2] = junction_copy_2

            graph.add_link(inputs[0], junction_copy_num_1)
            graph.add_link(inputs[1], junction_copy_num_2)

            if option_1_votes > option_2_votes:
                graph.add_link(junction_copy_num_1, outputs[0])
                graph.add_link(junction_copy_num_2, outputs[1])
            else:  # option 2
                graph.add_link(junction_copy_num_1, outputs[1])
                graph.add_link(junction_copy_num_2, outputs[0])
            table_row.append('applied')

        two_way_junctions_table.append(table_row)

    print_table(two_way_junctions_table, alignments='RLLRRRR', max_col_width=13,
                sub_colour={'too many alt': 'red', 'no support': 'red', 'too close': 'red',
                            'applied': 'green'})
    log.log('')


def simple_bridge_loops(graph, start_overlap_reads, end_overlap_reads, minimap_alignments,
                        read_dict):
    log.log('Simple loop bridges:')

    scoring_scheme = AlignmentScoringScheme('0,-1,-1,-1')

    loops = graph.find_all_simple_loops()
    for start, end, middle, repeat in loops:
        start_len = graph.segments[abs(start)].get_length()
        end_len = graph.segments[abs(end)].get_length()
        if start_len < settings.MIN_SEGMENT_LENGTH_FOR_SIMPLE_BRIDGING or \
                end_len < settings.MIN_SEGMENT_LENGTH_FOR_SIMPLE_BRIDGING:
            continue

        print()  # TEMP
        print()  # TEMP
        print()  # TEMP
        print()  # TEMP
        print(start, end, middle, repeat)  # TEMP

        forward_strand_reads = end_overlap_reads[start] & start_overlap_reads[end]
        reverse_strand_reads = end_overlap_reads[-end] & start_overlap_reads[-start]

        all_reads = list(forward_strand_reads) + list(reverse_strand_reads)
        strands = ['F'] * len(forward_strand_reads) + ['R'] * len(reverse_strand_reads)

        # This dictionary will collect the votes. The key is the number of times through the loop
        # and the value is the vote count. Votes for -1 times through the loop occur for reads
        # which don't conform to the loop assumption.
        votes = defaultdict(int)

        # We'll try a range of repeat counts, using the segment depth to guide how high we should
        # go.
        mean_start_end_depth = (graph.segments[abs(start)].depth +
                                graph.segments[abs(end)].depth) / 2
        best_repeat_guess = int(round(graph.segments[abs(middle)].depth / mean_start_end_depth))
        best_repeat_guess = max(1, best_repeat_guess)
        print('best_repeat_guess =', best_repeat_guess)  # TEMP
        max_tested_repeat = (best_repeat_guess + 1) * 2
        print('max_tested_repeat =', max_tested_repeat)  # TEMP

        for read, strand in zip(all_reads, strands):
            print(read, strand)  # TEMP

            if strand == 'F':
                s, e, m, r = start, end, middle, repeat
            else:  # strand == 'R'
                s, e, m, r = -end, -start, -middle, -repeat
            alignments = minimap_alignments[read]

            last_index_of_start = -1
            for i, a in enumerate(alignments):
                if a.get_signed_ref_name() == str(s):
                    last_index_of_start = i
            first_index_of_end = -1
            for i in range(last_index_of_start + 1, len(alignments)):
                a = alignments[i]
                if a.get_signed_ref_name() == str(e):
                    first_index_of_end = i
                    break

            # We should now have the indices of the alignments around the repeat.
            assert(last_index_of_start != -1)
            assert(first_index_of_end != -1)

            # If there are any alignment in between the start and end segments, they are only
            # allowed to be the middle or repeat segments.
            bad_middle = False
            for i in range(last_index_of_start+1, first_index_of_end):
                ref_name = alignments[i].get_signed_ref_name()
                if ref_name != str(m) and ref_name != str(r):
                    bad_middle = True
            if bad_middle:
                votes[-1] += 1  # vote for bad read
                continue

            start_alignment = alignments[last_index_of_start]
            end_alignment = alignments[first_index_of_end]
            print(start_alignment)  # TEMP
            print(end_alignment)  # TEMP

            # Now that we have the alignments, we can extract the relevant part of the read...
            read_start_pos, read_end_pos = start_alignment.read_start, end_alignment.read_end
            read_seq = read_dict[read].sequence[read_start_pos:read_end_pos]

            # ... and the relevant parts of the start/end segments.
            if start_alignment.read_strand == '+':
                start_seg_start_pos = start_alignment.ref_start
            else:  # start_alignment.read_strand == '-'
                start_seg_start_pos = start_alignment.ref_length - start_alignment.ref_end
            if end_alignment.read_strand == '+':
                end_seg_end_pos = end_alignment.ref_end
            else:  # end_alignment.read_strand == '-'
                end_seg_end_pos = end_alignment.ref_length - end_alignment.ref_start
            start_seg_seq = graph.seq_from_signed_seg_num(s)[start_seg_start_pos:]
            end_seg_seq = graph.seq_from_signed_seg_num(e)[:end_seg_end_pos]

            print('READ SEQ:  ', read_seq)  # TEMP

            middle_seq = graph.seq_from_signed_seg_num(m)
            repeat_seq = graph.seq_from_signed_seg_num(r)

            for repeat_count in range(0, max_tested_repeat + 1):
                test_seq = start_seg_seq + repeat_seq
                for _ in range(repeat_count):
                    test_seq += middle_seq + repeat_seq
                test_seq += end_seg_seq
                print('TEST SEQ ' + str(repeat_count) + ':', test_seq)  # TEMP
                alignment_result = fully_global_alignment(read_seq, test_seq, scoring_scheme, True,
                                                          settings.SIMPLE_REPEAT_BRIDGING_BAND_SIZE)
                if alignment_result:
                    seqan_parts = alignment_result.split(',', 9)
                    test_seq_score = int(seqan_parts[6])
                else:
                    test_seq_score = 0

                print('TEST SEQ ' + str(repeat_count) + ' SCORE:', test_seq_score)  # TEMP

            print()  # TEMP





