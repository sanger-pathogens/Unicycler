"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

Long read bridges are the big important type of bridge for a hybrid Unicycler assembly. They are
made using long reads which align to multiple segments in the graph.

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import time
import math
from .bridge_common import get_bridge_str, get_mean_depth, get_depth_agreement_factor
from .misc import float_to_str
from . import settings
from . import log
from .path_finding import get_best_paths_for_seq


class MiniasmBridge(object):
    """
    This class describes a bridge created from long read alignments.
    """
    def __init__(self, graph, start, end, bridge_sequence, start_overlap, end_overlap,
                 scoring_scheme):

        # The numbers of the two single copy segments which are being bridged.
        self.start_segment = start
        self.end_segment = end

        print('MINIASM BRIDGE: ' + str(start) + ' to ' + str(end))  # TEMP
        print('  BRIDGE SEQ LENGTH: ' + str(len(bridge_sequence)))  # TEMP

        # In some miniasm bridges where the contigs are very close, there will be overlap between
        # the contigs and the bridge. I.e. when the bridge is applied, contig sequence will be
        # trimmed off.
        self.start_overlap = start_overlap
        self.end_overlap = end_overlap
        print('  OVERLAPS: ' + str(start_overlap) + ', ' + str(end_overlap))  # TEMP

        # The bridge depth, a weighted mean of the start and end depths.
        self.depth = get_mean_depth(graph.segments[abs(self.start_segment)],
                                    graph.segments[abs(self.end_segment)], graph)

        # Look for a graph path corresponding to the bridge sequence.
        target_path_length = len(bridge_sequence)
        path_start_time = time.time()
        self.all_paths, progressive_path_search = \
            get_best_paths_for_seq(graph, self.start_segment, self.end_segment, target_path_length,
                                   bridge_sequence, scoring_scheme, 0.9)
        path_time = time.time() - path_start_time
        if self.all_paths:
            self.graph_path = self.all_paths[0][0]
            scaled_score = self.all_paths[0][3]
        else:
            self.graph_path = []
            scaled_score = 0.0
        print('  PATH SEARCH TIME: ' + str(path_time))  # TEMP
        print('  BEST PATH: ' + str(self.graph_path))  # TEMP
        print('  PATH SCALED SCORE: ' + str(scaled_score))  # TEMP

        # If a very good match was found, use the graph path and give the bridge a very high
        # quality score.
        if scaled_score > settings.MINIASM_BRIDGE_SCALED_SCORE_TO_USE_GRAPH_PATH:
            self.bridge_sequence = graph.get_path_sequence(self.graph_path)
            self.quality = settings.MINIASM_BRIDGE_QUAL_WITH_GRAPH_PATH
            print('  USING PATH')  # TEMP

        # Otherwise, just use the miniasm bridge sequence for the bridge.
        else:
            self.bridge_sequence = bridge_sequence
            if graph.ends_with_dead_end(self.start_segment) or \
                    graph.starts_with_dead_end(self.end_segment):
                self.quality = settings.MINIASM_BRIDGE_QUAL_WITH_DEAD_END
            else:
                self.quality = settings.MINIASM_BRIDGE_QUAL_WITHOUT_PATH_OR_DEAD_END
            print('  NOT USING PATH')  # TEMP

        # Depth agreement affects bridge quality.
        start_seg = graph.segments[abs(self.start_segment)]
        end_seg = graph.segments[abs(self.end_segment)]
        self.quality *= get_depth_agreement_factor(start_seg.depth, end_seg.depth)

        # Bridge length affects quality too: short bridges are better.
        bridge_len = max(0, len(self.bridge_sequence))
        half_qual_len = settings.MINIASM_BRIDGE_HALF_QUAL_LENGTH
        self.quality *= half_qual_len / (bridge_len + half_qual_len)

        self.quality = 100.0 * math.sqrt(self.quality)

        print('  FINAL QUALITY: ' + str(self.quality))  # TEMP
        print('')  # TEMP

        self.segments_reduced_depth = []

    def __repr__(self):
        return 'miniasm bridge: ' + get_bridge_str(self) + \
               ' (quality = ' + float_to_str(self.quality, 2) + ')'

    @staticmethod
    def get_type_score():
        """
        Returns a score indicating the relative importance of the bridge types:
        LongReadBridge = 2, SpadesContigBridge = 1, LoopUnrollingBridge = 0
        """
        return 2

    @staticmethod
    def get_type_name():
        """
        Returns the name of the bridge type.
        """
        return 'miniasm'


def create_miniasm_bridges(graph, string_graph, anchor_segments, scoring_scheme):
    """
    Makes bridges between single copy segments using the miniasm string graph.
    """
    log.log_section_header('Creating miniasm/Racon bridges')
    log.log_explanation('Now that the miniasm/Racon string graph is complete, Unicycler will '
                        'use it to build bridges between single-copy contigs.', verbosity=1)
    bridges = []
    seg_nums_to_bridge = set(x.number for x in anchor_segments)

    string_graph_bridge_segments = [x for x in string_graph.segments
                                    if x.startswith('BRIDGE_') or
                                    x.startswith('OVERLAPPING_BRIDGE_')]

    for bridge_seg_name in string_graph_bridge_segments:
        bridge_seg = string_graph.segments[bridge_seg_name]
        pos_seg_name = bridge_seg_name + '+'
        preceding_segments = string_graph.get_preceding_segments(pos_seg_name)
        following_segments = string_graph.get_following_segments(pos_seg_name)
        if len(preceding_segments) != 1:
            continue
        if len(following_segments) != 1:
            continue

        preceding_seg_name = preceding_segments[0]
        following_seg_name = following_segments[0]
        assert preceding_seg_name.startswith('CONTIG_')
        assert following_seg_name.startswith('CONTIG_')
        first_link = string_graph.links[(preceding_seg_name, pos_seg_name)]
        second_link = string_graph.links[(pos_seg_name, following_seg_name)]
        preceding_seg_name = preceding_seg_name[7:]
        following_seg_name = following_seg_name[7:]
        preceding_segment_number = int(preceding_seg_name[:-1]) * \
                                   (1 if preceding_seg_name[-1] == '+' else -1)
        following_segment_number = int(following_seg_name[:-1]) * \
                                   (1 if following_seg_name[-1] == '+' else -1)
        assert abs(preceding_segment_number) in seg_nums_to_bridge
        assert abs(following_segment_number) in seg_nums_to_bridge

        start_overlap = first_link.seg_1_overlap
        end_overlap = second_link.seg_2_overlap

        bridges.append(MiniasmBridge(graph, preceding_segment_number, following_segment_number,
                                     bridge_seg.forward_sequence, start_overlap, end_overlap,
                                     scoring_scheme))
    return bridges
