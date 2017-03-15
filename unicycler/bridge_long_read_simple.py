"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module takes care of simple long read bridging. This is where parts of the graph have a very
straightforward repeat and minimap alignments can conclusively support one over the alternatives.
Simple long read bridging is done before other, more complex bridging processes.

Specifically, there are two kinds of simple bridging:
 * two way junctions
 * simple loops

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
from multiprocessing.dummy import Pool as ThreadPool
from .cpp_wrappers import fully_global_alignment
from .minimap_alignment import align_long_reads_to_assembly_graph, build_start_end_overlap_sets
from .assembly_graph_segment import Segment
from .assembly_graph_copy_depth import determine_copy_depth
from .misc import print_table, get_right_arrow
from . import log
from . import settings


def apply_simple_long_read_bridges(graph, out_dir, keep, threads, read_dict, long_read_filename,
                                   scoring_scheme):
    """
    Look for and apply simple long read bridges to the graph.
    """
    log.log_section_header('Simple repeat bridging with minimap alignments')
    log.log_explanation('Unicycler uses long read alignments (from minimap) to resolve simple '
                        'repeat structures in the graph. This simplifies the graph before more '
                        'complex bridging approaches are applied in later steps.')

    bridging_dir = os.path.join(out_dir, 'simple_bridging')
    if not os.path.exists(bridging_dir):
        os.makedirs(bridging_dir)

    minimap_alignments = align_long_reads_to_assembly_graph(graph, long_read_filename,
                                                            bridging_dir, threads)
    start_overlap_reads, end_overlap_reads = build_start_end_overlap_sets(minimap_alignments)

    simple_bridge_two_way_junctions(graph, start_overlap_reads, end_overlap_reads,
                                    minimap_alignments)
    simple_bridge_loops(graph, start_overlap_reads, end_overlap_reads, minimap_alignments,
                        read_dict, scoring_scheme, threads)
    log.log('', 2)

    graph.merge_all_possible(None, 2)
    determine_copy_depth(graph)

    if keep < 3:
        shutil.rmtree(bridging_dir)
    return graph


def simple_bridge_two_way_junctions(graph, start_overlap_reads, end_overlap_reads,
                                    minimap_alignments):
    c_with_arrows = get_right_arrow() + 'C' + get_right_arrow()
    log.log_explanation('Two-way junctions are defined as cases where two graph contigs (A and B) '
                        'join together (C) and then split apart again (D and E). This usually '
                        'represents a simple 2-copy repeat, and there are two possible options for '
                        'its resolution: (A' + c_with_arrows + 'D and B' + c_with_arrows + 'E) or '
                        '(A' + c_with_arrows + 'E and B' + c_with_arrows + 'D). '
                        'Each read which spans such a junction gets to "vote" for option 1, '
                        'option 2 or neither. If the reads\' votes show option 1 or 2 as a clear '
                        'winner, Unicycler resolves the graph at that junction.')

    two_way_junctions_table = [['Junction', 'Option 1', 'Option 2', 'Op. 1 votes', 'Op. 2 votes',
                                'Neither votes', 'Result']]

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

        spaced_arrow = ' ' + get_right_arrow() + ' '
        option_1_1_str = spaced_arrow.join(str(x) for x in [inputs[0], junction, outputs[0]])
        option_1_2_str = spaced_arrow.join(str(x) for x in [inputs[1], junction, outputs[1]])
        table_row.append(option_1_1_str + ', ' + option_1_2_str)
        option_2_1_str = spaced_arrow.join(str(x) for x in [inputs[0], junction, outputs[1]])
        option_2_2_str = spaced_arrow.join(str(x) for x in [inputs[1], junction, outputs[0]])
        table_row.append(option_2_1_str + ', ' + option_2_2_str)

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

    max_op_1_len = max(max(len(y) for y in x[1].split(', ')) for x in two_way_junctions_table) + 1
    max_op_2_len = max(max(len(y) for y in x[2].split(', ')) for x in two_way_junctions_table) + 1
    max_result_len = max(len(x[-1]) for x in two_way_junctions_table)
    print_table(two_way_junctions_table, alignments='RLLRRRR', left_align_header=False, indent=0,
                fixed_col_widths=[8, max_op_1_len, max_op_2_len, 5, 5, 7, max_result_len],
                sub_colour={'too many alt': 'red', 'no support': 'red', 'too close': 'red',
                            'applied': 'green'})
    log.log('')


def simple_bridge_loops(graph, start_overlap_reads, end_overlap_reads, minimap_alignments,
                        read_dict, scoring_scheme, threads):
    ra = get_right_arrow()
    zero_loops = 'A' + ra + 'C' + ra + 'B'
    one_loop = 'A' + ra + 'C' + ra + 'D' + ra + 'C' + ra + 'B'
    two_loops = 'A' + ra + 'C' + ra + 'D' + ra + 'C' + ra + 'D' + ra + 'C' + ra + 'B'
    log.log_explanation('Simple loops are defined as parts of the graph where two contigs (A and '
                        'B) are connected via a repeat (C) which loops back to itself (via D). It '
                        'is possible to traverse the loop zero times (' + zero_loops + '), one '
                        'time (' + one_loop + '), two times (' + two_loops + '), etc. '
                        'Long reads which span the loop inform which is the correct number of '
                        'times through. In this step, such reads are found and each is aligned '
                        'against alternative loop counts. A reads casts its "vote" for the loop '
                        'count it agrees best with, and Unicycler resolves the graph using the '
                        'most voted for count.')

    col_widths = [5, 6, 6, 5, 5, 18, 17]
    loop_table_header = ['Start', 'Repeat', 'Middle', 'End', 'Read count', 'Read votes', 'Result']
    print_table([loop_table_header], fixed_col_widths=col_widths, left_align_header=False,
                alignments='RRRRRLL', indent=0)

    loops = sorted(graph.find_all_simple_loops())
    for start, end, middle, repeat in loops:
        start_len = graph.segments[abs(start)].get_length()
        end_len = graph.segments[abs(end)].get_length()
        if start_len < settings.MIN_SEGMENT_LENGTH_FOR_SIMPLE_BRIDGING or \
                end_len < settings.MIN_SEGMENT_LENGTH_FOR_SIMPLE_BRIDGING:
            continue

        loop_table_row = [start, repeat, middle, end]

        forward_strand_reads = end_overlap_reads[start] & start_overlap_reads[end]
        reverse_strand_reads = end_overlap_reads[-end] & start_overlap_reads[-start]

        all_reads = list(forward_strand_reads) + list(reverse_strand_reads)
        strands = ['F'] * len(forward_strand_reads) + ['R'] * len(reverse_strand_reads)
        loop_table_row.append(len(all_reads))

        # This dictionary will collect the votes. The key is the number of times through the loop
        # and the value is the vote count. Votes for -1 times through the loop occur for reads
        # which don't conform to the loop assumption.
        votes = defaultdict(int)

        # We'll try a range of repeat counts. The segment depth gives us a first guess as to the
        # repeat count, which guides how high we should test.
        mean_start_end_depth = (graph.segments[abs(start)].depth +
                                graph.segments[abs(end)].depth) / 2
        best_repeat_guess = int(round(graph.segments[abs(middle)].depth / mean_start_end_depth))
        best_repeat_guess = max(1, best_repeat_guess)
        max_tested_loop_count = (best_repeat_guess + 1) * 2

        # Use a simple loop if we only have one thread.
        if threads == 1:
            for read, strand in zip(all_reads, strands):
                vote = get_read_loop_vote(start, end, middle, repeat, strand, minimap_alignments,
                                          read, read_dict, graph, max_tested_loop_count,
                                          scoring_scheme)
                votes[vote] += 1

        # Use a thread pool if we have more than one thread.
        else:
            pool = ThreadPool(threads)
            arg_list = []
            for read, strand in zip(all_reads, strands):
                arg_list.append((start, end, middle, repeat, strand, minimap_alignments,
                                 read, read_dict, graph, max_tested_loop_count, scoring_scheme))
            for vote in pool.imap_unordered(get_read_loop_vote_one_arg, arg_list):
                votes[vote] += 1

        # Format the vote totals nicely for the table.
        vote_str = ''
        for loop_count in sorted(votes.keys()):
            if loop_count == -1:
                vote_str += 'bad: '
            elif loop_count == 1:
                vote_str += '1 loop: '
            else:
                vote_str += str(loop_count) + ' loops: '
            vote_count = votes[loop_count]
            vote_str += str(vote_count) + ' vote' + ('s' if vote_count != 1 else '') + '    '
        loop_table_row.append(vote_str.strip())

        # Determine the repeat count which wins!
        applied_message = None
        results = sorted(list(votes.items()), key=lambda x: x[1], reverse=True)
        if not results:
            loop_table_row.append('no reads')
        else:
            winning_loop_count = results[0][0]
            if len(results) == 1:
                second_best_ratio = 0.0
            else:
                second_best_ratio = results[1][1] / results[0][1]
            if winning_loop_count == -1:
                loop_table_row.append('too many bad')
            elif second_best_ratio > 0.5:
                loop_table_row.append('too close')
            else:
                applied_message = 'applied: ' + str(winning_loop_count) + ' loop' + \
                                  ('' if winning_loop_count == 1 else 's')
                loop_table_row.append(applied_message)

                # If we got here, then we're good to bridge!
                middle_seq = graph.seq_from_signed_seg_num(middle)
                repeat_seq = graph.seq_from_signed_seg_num(repeat)
                new_seg_seq = repeat_seq
                for _ in range(winning_loop_count):
                    new_seg_seq += middle_seq + repeat_seq
                new_seg_num = graph.get_next_available_seg_number()
                new_seg = Segment(new_seg_num, mean_start_end_depth, new_seg_seq, True)
                graph.segments[new_seg_num] = new_seg
                graph.add_link(start, new_seg_num)
                graph.add_link(new_seg_num, end)

                # If we've traversed the loop, we can delete the segments. If we've haven't, then
                # we leave them as a separate circular sequence.
                if winning_loop_count > 0:
                    graph.remove_segments([abs(middle), abs(repeat)])
                else:
                    graph.remove_link(start, repeat)
                    graph.remove_link(repeat, end)

        sub_colours = {'too many bad': 'red', 'no reads': 'red', 'too close': 'red'}
        if applied_message:
            sub_colours[applied_message] = 'green'
        print_table([loop_table_row], fixed_col_widths=col_widths, header_format='normal',
                    alignments='RRRRRLL', left_align_header=False, bottom_align_header=False,
                    sub_colour=sub_colours, indent=0)


def get_read_loop_vote_one_arg(all_args):
    start, end, middle, repeat, strand, minimap_alignments, read, read_dict, graph, \
        max_tested_loop_count, scoring_scheme = all_args
    return get_read_loop_vote(start, end, middle, repeat, strand, minimap_alignments, read, read_dict,
                       graph, max_tested_loop_count, scoring_scheme)


def get_read_loop_vote(start, end, middle, repeat, strand, minimap_alignments, read, read_dict,
                       graph, max_tested_loop_count, scoring_scheme):
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
    assert (last_index_of_start != -1)
    assert (first_index_of_end != -1)

    # If there are any alignments in between the start and end segments, they are only
    # allowed to be the middle or repeat segments.
    bad_middle = False
    for i in range(last_index_of_start + 1, first_index_of_end):
        ref_name = alignments[i].get_signed_ref_name()
        if ref_name != str(m) and ref_name != str(r):
            bad_middle = True
    if bad_middle:
        return -1  # vote for bad read

    start_alignment = alignments[last_index_of_start]
    end_alignment = alignments[first_index_of_end]

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

    middle_seq = graph.seq_from_signed_seg_num(m)
    repeat_seq = graph.seq_from_signed_seg_num(r)

    best_score = None
    best_count = None
    for loop_count in range(0, max_tested_loop_count + 1):
        test_seq = start_seg_seq + repeat_seq
        for _ in range(loop_count):
            test_seq += middle_seq + repeat_seq
        test_seq += end_seg_seq
        alignment_result = fully_global_alignment(read_seq, test_seq, scoring_scheme, True,
                                                  settings.SIMPLE_REPEAT_BRIDGING_BAND_SIZE)
        if alignment_result:
            seqan_parts = alignment_result.split(',', 9)
            test_seq_score = int(seqan_parts[6])
            if best_score is None or test_seq_score > best_score:
                best_score = test_seq_score
                best_count = loop_count

    # This read now casts its vote for the best repeat count!
    if best_count is not None:
        return best_count
