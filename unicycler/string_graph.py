"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This module describes a miniasm string graph and related functions. Unlike the AssemblyGraph class
(which hold a SPAdes-generated graph), a miniasm string graph does not have constant overlaps.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import statistics
from collections import defaultdict
from .misc import reverse_complement, get_num_agreement
from .assembly_graph import build_reverse_links
from . import settings
from . import log


class StringGraph(object):

    def __init__(self, filename, mean_read_quals):
        self.segments = {}                      # unsigned seg name -> StringGraphSegment
        self.forward_links = defaultdict(list)  # signed seg name -> list of signed segment name
        self.reverse_links = defaultdict(list)  # signed seg name <- list of signed segment name
        self.links = {}                         # tuple (start, end) -> StringGraphLink

        # Load in the segments.
        with open(filename, 'rt') as gfa_file:
            for line in gfa_file:
                if line.startswith('S'):
                    line_parts = line.strip().split('\t')
                    name = line_parts[1]
                    sequence = line_parts[2]
                    self.segments[name] = StringGraphSegment(name, sequence, mean_read_quals)

        # Load in the links.
        with open(filename, 'rt') as gfa_file:
            for line in gfa_file:
                if line.startswith('L'):
                    line_parts = line.strip().split('\t')
                    signed_name_1 = line_parts[1] + line_parts[2]
                    signed_name_2 = line_parts[3] + line_parts[4]
                    self.forward_links[signed_name_1].append(signed_name_2)

                    link_tuple = (signed_name_1, signed_name_2)
                    if link_tuple not in self.links:
                        self.links[link_tuple] = StringGraphLink(signed_name_1, signed_name_2)
                    seg_1_to_seg_2_overlap = int(line_parts[5][:-1])
                    self.links[link_tuple].seg_1_overlap = seg_1_to_seg_2_overlap

                    rev_name_1 = flip_segment_name(signed_name_1)
                    rev_name_2 = flip_segment_name(signed_name_2)
                    rev_link_tuple = (rev_name_2, rev_name_1)
                    if rev_link_tuple not in self.links:
                        self.links[rev_link_tuple] = StringGraphLink(rev_name_2, rev_name_1)
                    self.links[rev_link_tuple].seg_2_overlap = seg_1_to_seg_2_overlap
            self.reverse_links = build_reverse_links(self.forward_links)

    def save_to_gfa(self, filename, verbosity=1, newline=False):
        """
        Saves whole graph to a GFA file.
        """
        log.log(('\n' if newline else '') + 'Saving ' + filename, verbosity)
        with open(filename, 'w') as gfa:
            for segment in sorted(self.segments.values(), key=lambda x: x.full_name):
                gfa.write(segment.gfa_segment_line())
            for link in sorted(self.links.keys()):
                gfa.write(self.links[link].gfa_link_line())

    def get_preceding_segments(self, seg_name):
        if seg_name not in self.reverse_links:
            return []
        return self.reverse_links[seg_name]

    def get_following_segments(self, seg_name):
        if seg_name not in self.forward_links:
            return []
        return self.forward_links[seg_name]

    def add_link(self, start, end, overlap_1, overlap_2):
        """
        Adds a link to the graph in all necessary ways: forward and reverse, and for reverse
        complements too.
        """
        rev_start = flip_segment_name(start)
        rev_end = flip_segment_name(end)

        if start not in self.forward_links:
            self.forward_links[start] = []
        if end not in self.forward_links[start]:
            self.forward_links[start].append(end)

        if end not in self.reverse_links:
            self.reverse_links[end] = []
        if start not in self.reverse_links[end]:
            self.reverse_links[end].append(start)

        if rev_start not in self.reverse_links:
            self.reverse_links[rev_start] = []
        if rev_end not in self.reverse_links[rev_start]:
            self.reverse_links[rev_start].append(rev_end)

        if rev_end not in self.forward_links:
            self.forward_links[rev_end] = []
        if rev_start not in self.forward_links[rev_end]:
            self.forward_links[rev_end].append(rev_start)

        link_tuple = (start, end)
        self.links[link_tuple] = StringGraphLink(start, end)
        self.links[link_tuple].seg_1_overlap = overlap_1
        self.links[link_tuple].seg_2_overlap = overlap_2

        rev_link_tuple = (rev_end, rev_start)
        self.links[rev_link_tuple] = StringGraphLink(rev_end, rev_start)
        self.links[rev_link_tuple].seg_1_overlap = overlap_2
        self.links[rev_link_tuple].seg_2_overlap = overlap_1

    def remove_segment(self, seg_name_to_remove):
        """
        Removes a segment from the graph and all of its related links.
        """
        def remove_signed_segment(graph, seg_name):
            for preceding_seg_name in graph.get_preceding_segments(seg_name):
                del graph.links[(preceding_seg_name, seg_name)]
                graph.forward_links[preceding_seg_name].remove(seg_name)
            for following_seg_name in graph.get_following_segments(seg_name):
                del graph.links[(seg_name, following_seg_name)]
                graph.reverse_links[following_seg_name].remove(seg_name)
            graph.forward_links.pop(seg_name, None)
            graph.reverse_links.pop(seg_name, None)

        remove_signed_segment(self, seg_name_to_remove + '+')
        remove_signed_segment(self, seg_name_to_remove + '-')
        self.segments.pop(seg_name_to_remove, None)

    def remove_non_bridging_paths(self):
        """
        Removes any segments from the graph which are not on a simple path between two contigs, as
        these are not useful for bridging.
        """
        segments_to_remove = []
        for seg_name, segment in self.segments.items():
            if not self.segment_leads_directly_to_contig_in_both_directions(seg_name):
                segments_to_remove.append(seg_name)
        for seg_name in segments_to_remove:
            self.remove_segment(seg_name)

    def segment_leads_directly_to_contig_in_both_directions(self, seg_name):
        if self.segments[seg_name].contig:
            return True
        return (self.segment_leads_directly_to_contig(seg_name + '+') and
                self.segment_leads_directly_to_contig(seg_name + '-'))

    def segment_leads_directly_to_contig(self, signed_seg_name):
        """
        Tests whether a given segment leads to a contig via a simple unbranching path. Only tests
        in a single direction.
        """
        current_seg_name = signed_seg_name
        while True:
            following_segments = self.get_following_segments(current_seg_name)
            preceding_segments = self.get_preceding_segments(current_seg_name)
            if len(following_segments) != 1 or len(preceding_segments) != 1:
                return False
            if self.segments[get_unsigned_seg_name(current_seg_name)].contig:
                return True
            current_seg_name = following_segments[0]

    def get_bridging_paths(self):
        """
        Returns a list of all bridging paths. The contigs being bridged are included at the start
        and end of each path.
        """
        paths = []
        used_segments = set()
        for seg_name, segment in self.segments.items():
            if not segment.contig and seg_name not in used_segments and \
                    self.segment_leads_directly_to_contig_in_both_directions(seg_name):
                starting_seg = seg_name + '+'
                current_seg = starting_seg
                path = [current_seg]
                while True:
                    current_seg = self.get_following_segments(current_seg)[0]
                    path.append(current_seg)
                    if self.segments[get_unsigned_seg_name(current_seg)].contig:
                        break
                current_seg = starting_seg
                while True:
                    current_seg = self.get_preceding_segments(current_seg)[0]
                    path.insert(0, current_seg)
                    if self.segments[get_unsigned_seg_name(current_seg)].contig:
                        break
                for seg in path:
                    used_segments.add(get_unsigned_seg_name(seg))
                paths.append(path)
        return paths

    def simplify_bridges(self, before_transitive_reduction):
        """
        This function takes cases where there is a simple path from one contig to the next, and it
        replaces it with a single spanning read.
        """
        paths = self.get_bridging_paths()
        for path in paths:
            print('')  # TEMP
            print('SIMPLIFY BRIDGE:', path)  # TEMP
            assert len(path) >= 3
            contig_1 = path[0]
            contig_2 = path[-1]
            middle = path[1:-1]

            # If there is only one read joining the contigs, then there's nothing to simplify.
            if len(middle) == 1:
                continue

            # If there are multiple reads, we hope to find one that connects the two contigs.
            # Search them in order of highest to lowest quality, seeing if the links we need are
            # in the pre-transitive reduction graph.
            single_bridge_read = None
            reads_by_qual = sorted([x for x in middle], reverse=True,
                               key=lambda x: self.segments[get_unsigned_seg_name(x)].qual)
            for read in reads_by_qual:  # TEMP
                print(read + ':', self.segments[get_unsigned_seg_name(read)].qual)  # TEMP
            for read in reads_by_qual:
                if (contig_1, read) in before_transitive_reduction.links and \
                        (read, contig_2) in before_transitive_reduction.links:
                    single_bridge_read = read
                    break

            print('SINGLE BRIDGE READ:', single_bridge_read)  # TEMP
            if single_bridge_read is not None:

                # Delete all of the other segments in the bridge.
                for read in [x for x in middle if x != single_bridge_read]:
                    self.remove_segment(get_unsigned_seg_name(read))

                # Create links between the surviving segment and the contigs.
                link_1_tuple = (contig_1, single_bridge_read)
                link_2_tuple = (single_bridge_read, contig_2)
                if link_1_tuple not in self.links:
                    link_1 = before_transitive_reduction.links[link_1_tuple]
                    self.add_link(link_1_tuple[0], link_1_tuple[1],
                                  link_1.seg_1_overlap, link_1.seg_2_overlap)
                if link_2_tuple not in self.links:
                    link_2 = before_transitive_reduction.links[link_2_tuple]
                    self.add_link(link_2_tuple[0], link_2_tuple[1],
                                  link_2.seg_1_overlap, link_2.seg_2_overlap)


    def remove_overlaps(self, before_transitive_reduction):
        """
        Removes reads and trims reads to get rid of graph overlap:
          * Selectively removes overlaps from lower quality sequences first.
          * Keeps as much contig sequence as possible.

        before_transitive_reduction is the string graph at an earlier state where there are still
        transitive links. This helps us to find the correct overlaps when removing segments.
        This function assumes that remove_non_bridging_paths has been run and we only have simple,
        unbranching paths to deal with.
        """
        segments_by_quality = sorted([x for x in self.segments],
                                     key=lambda x: self.segments[x].qual)
        print('\n\n\n')  # TEMP
        for seg_name in segments_by_quality:
            seg = self.segments[seg_name]

            # We don't trim or remove contigs.
            if seg.contig:
                continue

            seg_len = seg.get_length()
            pos_seg_name = seg_name + '+'
            neg_seg_name = seg_name + '-'

            print('')  # TEMP
            print('REMOVE OVERLAPS:', pos_seg_name)  # TEMP

            preceding_segments = self.get_preceding_segments(pos_seg_name)
            following_segments = self.get_following_segments(pos_seg_name)

            print('  PRECEDING SEGMENTS:', preceding_segments)
            print('  FOLLOWING SEGMENTS:', following_segments)

            if len(preceding_segments) == 1 and len(following_segments) == 1:
                preceding_seg_name = preceding_segments[0]
                following_seg_name = following_segments[0]
                start_link = self.links[(preceding_seg_name, pos_seg_name)]
                rev_start_link = self.links[(neg_seg_name, flip_segment_name(preceding_seg_name))]
                end_link = self.links[(pos_seg_name, following_seg_name)]
                rev_end_link = self.links[(flip_segment_name(following_seg_name), neg_seg_name)]

                start_overlap = start_link.seg_2_overlap
                end_overlap = end_link.seg_1_overlap
                print('start_overlap:', start_overlap)  # TEMP
                print('end_overlap:  ', end_overlap)  # TEMP

                # If there aren't any overlaps, then we've nothing to do!
                if start_overlap == 0 and end_overlap == 0:
                    print('NO OVERLAPS')  # TEMP
                    continue

                # If the start and end overlap sum to less than the length of the segment, then we
                # trim off the overlaps and leave the segment in the middle.
                if start_overlap + end_overlap < seg_len:
                    print('TRIM ENDS')  # TEMP
                    print('  LENGTH BEFORE:', len(seg.forward_sequence))  # TEMP
                    if end_overlap > 0:
                        seg.forward_sequence = seg.forward_sequence[start_overlap:-end_overlap]
                    else:
                        seg.forward_sequence = seg.forward_sequence[start_overlap:]
                    seg.reverse_sequence = reverse_complement(seg.forward_sequence)
                    start_link.seg_1_overlap = 0
                    start_link.seg_2_overlap = 0
                    rev_start_link.seg_1_overlap = 0
                    rev_start_link.seg_2_overlap = 0
                    end_link.seg_1_overlap = 0
                    end_link.seg_2_overlap = 0
                    rev_end_link.seg_1_overlap = 0
                    rev_end_link.seg_2_overlap = 0
                    print('  LENGTH AFTER:', len(seg.forward_sequence))  # TEMP

                # If the start and end overlap are more than the length of the segment, then we can
                # remove the segment entirely (because the preceding and following segments will
                # still overlap).
                else:
                    print('REMOVE WHOLE SEGMENT')  # TEMP
                    print('  NEW LINK: ' + preceding_seg_name + ' -> ' + following_seg_name)  # TEMP

                    start_overlap_ratio = start_link.seg_1_overlap / start_link.seg_2_overlap
                    end_overlap_ratio = end_link.seg_2_overlap / end_link.seg_1_overlap

                    # Make a guess for the new link overlaps.
                    new_link_overlap = -(seg_len - start_overlap - end_overlap)
                    overlap_1 = int(round(new_link_overlap * start_overlap_ratio))
                    overlap_2 = int(round(new_link_overlap * end_overlap_ratio))
                    print('  GUESS OVERLAP 1: ', overlap_1)  # TEMP
                    print('  GUESS OVERLAP 2: ', overlap_2)  # TEMP

                    # Look for the link in the before transitive reduction graph to find the
                    # exact overlaps, they are available.
                    link_tuple = (preceding_seg_name, following_seg_name)
                    if link_tuple in before_transitive_reduction.links:
                        exact_link = before_transitive_reduction.links[link_tuple]
                        overlap_1 = exact_link.seg_1_overlap
                        overlap_2 = exact_link.seg_2_overlap
                        print('  EXACT OVERLAP 1: ', exact_link.seg_1_overlap)  # TEMP
                        print('  EXACT OVERLAP 2: ', exact_link.seg_2_overlap)  # TEMP

                    # If we failed to find the link in the graph, we can go directly to the
                    # overlaps file and look there.
                    else:
                        # TO DO
                        # TO DO
                        # TO DO
                        # TO DO
                        a = 5  # TEMP
                        a += 1  # TEMP

                    self.add_link(preceding_seg_name, following_seg_name, overlap_1, overlap_2)
                    self.remove_segment(seg_name)

            print('')  # TEMP

        # We should now hopefully have an overlap-free graph!
        try:
            for seg_name in self.segments.keys():
                pos_seg_name = seg_name + '+'
                neg_seg_name = seg_name + '-'
                preceding_segments = self.get_preceding_segments(pos_seg_name)
                following_segments = self.get_following_segments(pos_seg_name)
                assert len(preceding_segments) < 2
                assert len(following_segments) < 2
                if len(preceding_segments) == 1:
                    preceding_seg_name = preceding_segments[0]
                    start_link = self.links[(preceding_seg_name, pos_seg_name)]
                    rev_start_link = self.links[(neg_seg_name, flip_segment_name(preceding_seg_name))]
                    assert start_link.seg_1_overlap == 0
                    assert start_link.seg_2_overlap == 0
                    assert rev_start_link.seg_1_overlap == 0
                    assert rev_start_link.seg_2_overlap == 0
                if len(following_segments) == 1:
                    following_seg_name = following_segments[0]
                    end_link = self.links[(pos_seg_name, following_seg_name)]
                    rev_end_link = self.links[(flip_segment_name(following_seg_name), neg_seg_name)]
                    assert end_link.seg_1_overlap == 0
                    assert end_link.seg_2_overlap == 0
                    assert rev_end_link.seg_1_overlap == 0
                    assert rev_end_link.seg_2_overlap == 0
            print('GRAPH IS OVERLAP FREE!!!!!')  # TEMP
        except AssertionError:
            print('GRAPH STILL HAS OVERLAPS - BOOOOOOOOOO')  # TEMP

    # def merge_reads(self):
    #     """
    #     This function takes any reads in a simple path and merges them together. It assumes that
    #     the graph is now overlap-free.
    #     """
    #     segments_in_paths = set()
    #     paths_to_merge = []
    #     for seg_name in self.segments.keys():
    #         if seg_name in segments_in_paths:
    #             continue
    #         path = self.get_simple_read_path(seg_name + '+')
    #         if path:
    #             segments_in_paths.update(get_unsigned_seg_name(x) for x in path)
    #         if len(path) > 1:
    #             paths_to_merge.append(path)
    #
    #     print('\n\nPATHS TO MERGE:')  # TEMP
    #     for path in paths_to_merge:
    #         print(path)  # TEMP
    #         merged_seg_name = 'MERGED_' + '_'.join(get_unsigned_seg_name(x) for x in path)
    #         merged_seg_seq = ''
    #         merged_seq_quals = []
    #         for path_seg in path:
    #             seq = self.seq_from_signed_seg_name(path_seg)
    #             merged_seq_quals.append(self.segments[get_unsigned_seg_name(path_seg)].qual)
    #             merged_seg_seq += seq
    #         merged_seq_qual = statistics.mean(merged_seq_quals)
    #         print('')  # TEMP
    #         self.segments[merged_seg_name] = StringGraphSegment(merged_seg_name, merged_seg_seq)
    #         self.segments[merged_seg_name].qual = merged_seq_qual
    #         pos_merged_seg_name = merged_seg_name + '+'
    #
    #         preceding_segments = self.get_preceding_segments(path[0])
    #         following_segments = self.get_following_segments(path[-1])
    #         assert len(preceding_segments) == 1
    #         assert len(following_segments) == 1
    #         preceding_segment = preceding_segments[0]
    #         following_segment = following_segments[0]
    #         assert preceding_segment.startswith('CONTIG_')
    #         assert following_segment.startswith('CONTIG_')
    #
    #         self.add_link(preceding_segment, pos_merged_seg_name, 0, 0)
    #         self.add_link(pos_merged_seg_name, following_segment, 0, 0)
    #
    #         for path_seg in path:
    #             self.remove_segment(get_unsigned_seg_name(path_seg))
    #
    # def get_simple_read_path(self, starting_seg):
    #     """
    #     Starting from the given segment, this simple paths consisting only of reads (not contigs).
    #     """
    #     if self.segments[get_unsigned_seg_name(starting_seg)].contig:
    #         return []
    #     simple_path = [starting_seg]
    #
    #     # Extend path forward as much as possible.
    #     current_seg_name = starting_seg
    #     while True:
    #         following_segments = self.get_following_segments(current_seg_name)
    #         if len(following_segments) != 1:
    #             break
    #         current_seg_name = following_segments[0]
    #         if self.segments[get_unsigned_seg_name(current_seg_name)].contig:
    #             break
    #         simple_path.append(current_seg_name)
    #
    #     # Extend path backward as much as possible.
    #     current_seg_name = starting_seg
    #     while True:
    #         preceding_segments = self.get_preceding_segments(current_seg_name)
    #         if len(preceding_segments) != 1:
    #             break
    #         current_seg_name = preceding_segments[0]
    #         if self.segments[get_unsigned_seg_name(current_seg_name)].contig:
    #             break
    #         simple_path = [current_seg_name] + simple_path
    #
    #     return simple_path
    #
    # def seq_from_signed_seg_name(self, signed_name):
    #     assert(signed_name.endswith('+') or signed_name.endswith('-'))
    #     unsigned_seg_name = get_unsigned_seg_name(signed_name)
    #     if signed_name.endswith('+'):
    #         return self.segments[unsigned_seg_name].forward_sequence
    #     else:
    #         return self.segments[unsigned_seg_name].reverse_sequence


class StringGraphSegment(object):

    def __init__(self, full_name, sequence, mean_read_quals=None):
        self.full_name = full_name  # Has range at the end
        self.forward_sequence = sequence
        self.reverse_sequence = reverse_complement(sequence)

        # Miniasm trims reads and puts the start/end positions in the name.
        name_parts = full_name.rsplit(':', 1)
        assert(len(name_parts) == 2)
        self.short_name = name_parts[0]
        range = name_parts[1][:-1]
        self.start_pos, self.end_pos = (int(x) for x in range.split('-'))

        if self.short_name.startswith('CONTIG_'):
            self.contig = True
            self.qual = settings.CONTIG_READ_QSCORE
        elif mean_read_quals is not None:
            assert(self.short_name in mean_read_quals)  # if not a contig, it should be a real long read
            self.contig = False
            self.qual = mean_read_quals[self.short_name]

    def __repr__(self):
        if len(self.forward_sequence) > 6:
            seq_string = (self.forward_sequence[:3] + '...' + self.forward_sequence[-3:] + ', ' +
                          str(len(self.forward_sequence)) + ' bp')
        else:
            seq_string = self.forward_sequence
        return self.full_name + ' (' + seq_string + '), mean score = ' + str(self.qual)

    def get_length(self):
        return len(self.forward_sequence)

    def gfa_segment_line(self):
        return '\t'.join(['S', self.full_name, self.forward_sequence]) + '\n'


class StringGraphLink(object):

    def __init__(self, seg_1_signed_name, seg_2_signed_name):
        self.seg_1_signed_name = seg_1_signed_name
        self.seg_2_signed_name = seg_2_signed_name
        self.seg_1_overlap = None
        self.seg_2_overlap = None

    def __repr__(self):
        return (self.seg_1_signed_name + ' -> ' + self.seg_2_signed_name + ' (' +
                str(self.seg_1_overlap) + ', ' + str(self.seg_2_overlap) + ')')

    def gfa_link_line(self):
        seg_1_name = get_unsigned_seg_name(self.seg_1_signed_name)
        seg_1_sign = self.seg_1_signed_name[-1]
        seg_2_name = get_unsigned_seg_name(self.seg_2_signed_name)
        seg_2_sign = self.seg_2_signed_name[-1]
        overlap = str(self.seg_1_overlap) + 'M'
        return '\t'.join(['L', seg_1_name, seg_1_sign, seg_2_name, seg_2_sign, overlap]) + '\n'


def flip_segment_name(seg_name):
    assert(seg_name.endswith('+') or seg_name.endswith('-'))
    if seg_name.endswith('+'):
        return get_unsigned_seg_name(seg_name) + '-'
    else:
        return get_unsigned_seg_name(seg_name) + '+'

def get_unsigned_seg_name(seg_name):
    assert(seg_name.endswith('+') or seg_name.endswith('-'))
    return seg_name[:-1]
