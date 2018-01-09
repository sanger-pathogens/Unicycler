"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Unicycler

This file is part of Unicycler. Unicycler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Unicycler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Unicycler. If
not, see <http://www.gnu.org/licenses/>.
"""

import unittest
import os
import shutil
import unicycler.miniasm_assembly
import unicycler.log
import unicycler.assembly_graph
import unicycler.string_graph
import unicycler.alignment


def sequences_match_some_rotation(seq_1, seq_2):
    if len(seq_1) != len(seq_2):
        return False
    for i in range(len(seq_1)):
        rotated_seq_1 = seq_1[i:] + seq_1[:i]
        if rotated_seq_1 == seq_2:
            return True
    return False


def get_merged_string_graph_seqs(string_graph):
    unitig_graph = \
        unicycler.string_graph.merge_string_graph_segments_into_unitig_graph(string_graph, dict())
    return sorted([x.forward_sequence for x in unitig_graph.segments.values()],
                  key=lambda s: len(s), reverse=True)


class TestContigPlacement(unittest.TestCase):
    """
    These tests look at contig placement into a unitig graph. Specifically, it shouldn't matter
    what the rotation of the unitig is, you should get the same contig-placed graph regardless.
    """

    def setUp(self):
        self.working_dir = self.blast_dir = 'TEMP_' + str(os.getpid())
        if not os.path.exists(self.working_dir):
            os.makedirs(self.working_dir)

        assembly_graph_filename = os.path.join(os.path.dirname(__file__),
                                               'test_contig_placement_assembly_graph.gfa')
        self.assembly_graph = unicycler.assembly_graph.AssemblyGraph(assembly_graph_filename,
                                                                     None)
        self.seg_nums_to_bridge = {122, 124, 125, 126, 237, 239}
        self.scoring_scheme = unicycler.alignment.AlignmentScoringScheme('3,-6,-5,-2')
        unicycler.log.logger = unicycler.log.Log(log_filename=None, stdout_verbosity_level=0)

    def tearDown(self):
        if os.path.exists(self.working_dir):
            shutil.rmtree(self.working_dir)

    def test_contig_placement_1(self):
        """
        This covers the case where the contig sequences are in the middle of the unitig sequences.
        """
        unitig_graph_filename = os.path.join(os.path.dirname(__file__),
                                             'test_contig_placement_unitig_graph_1.gfa')
        unitig_graph = unicycler.string_graph.StringGraph(unitig_graph_filename)
        new_graph = unicycler.miniasm_assembly.place_contigs(self.working_dir, self.assembly_graph,
                                                             unitig_graph, 1, self.scoring_scheme,
                                                             self.seg_nums_to_bridge)
        self.assertTrue(len(new_graph.segments), 12)
        self.assertEqual(new_graph.get_total_segment_length(), 18391)

        merged_seqs = get_merged_string_graph_seqs(new_graph)
        self.assertTrue(sequences_match_some_rotation(merged_seqs[0],
                                                      unitig_graph.segments['1'].forward_sequence))
        self.assertTrue(sequences_match_some_rotation(merged_seqs[1],
                                                      unitig_graph.segments['2'].forward_sequence))
        self.assertTrue(sequences_match_some_rotation(merged_seqs[2],
                                                      unitig_graph.segments['3'].forward_sequence))
        self.assertTrue(sequences_match_some_rotation(merged_seqs[3],
                                                      unitig_graph.segments['4'].forward_sequence))


    def test_contig_placement_2(self):
        """
        This covers the case where a contig sequence spans the circular junction of the unitig
        sequences.
        """
        unitig_graph_filename = os.path.join(os.path.dirname(__file__),
                                             'test_contig_placement_unitig_graph_2.gfa')
        unitig_graph = unicycler.string_graph.StringGraph(unitig_graph_filename)
        new_graph = unicycler.miniasm_assembly.place_contigs(self.working_dir, self.assembly_graph,
                                                             unitig_graph, 1, self.scoring_scheme,
                                                             self.seg_nums_to_bridge)
        self.assertTrue(len(new_graph.segments), 12)
        self.assertEqual(new_graph.get_total_segment_length(), 18391)

        merged_seqs = get_merged_string_graph_seqs(new_graph)
        self.assertTrue(sequences_match_some_rotation(merged_seqs[0],
                                                      unitig_graph.segments['1'].forward_sequence))
        self.assertTrue(sequences_match_some_rotation(merged_seqs[1],
                                                      unitig_graph.segments['2'].forward_sequence))
        self.assertTrue(sequences_match_some_rotation(merged_seqs[2],
                                                      unitig_graph.segments['3'].forward_sequence))
        self.assertTrue(sequences_match_some_rotation(merged_seqs[3],
                                                      unitig_graph.segments['4'].forward_sequence))

    def test_contig_placement_3(self):
        """
        This covers the case where the a contig sequence starts right at the beginning of the
        unitig sequence.
        """
        unitig_graph_filename = os.path.join(os.path.dirname(__file__),
                                             'test_contig_placement_unitig_graph_3.gfa')
        unitig_graph = unicycler.string_graph.StringGraph(unitig_graph_filename)
        new_graph = unicycler.miniasm_assembly.place_contigs(self.working_dir, self.assembly_graph,
                                                             unitig_graph, 1, self.scoring_scheme,
                                                             self.seg_nums_to_bridge)
        self.assertTrue(len(new_graph.segments), 12)
        self.assertEqual(new_graph.get_total_segment_length(), 18391)

        merged_seqs = get_merged_string_graph_seqs(new_graph)
        self.assertTrue(sequences_match_some_rotation(merged_seqs[0],
                                                      unitig_graph.segments['1'].forward_sequence))
        self.assertTrue(sequences_match_some_rotation(merged_seqs[1],
                                                      unitig_graph.segments['2'].forward_sequence))
        self.assertTrue(sequences_match_some_rotation(merged_seqs[2],
                                                      unitig_graph.segments['3'].forward_sequence))
        self.assertTrue(sequences_match_some_rotation(merged_seqs[3],
                                                      unitig_graph.segments['4'].forward_sequence))

    def test_contig_placement_4(self):
        """
        This covers the case where the a contig sequence ends right at the end of the unitig
        sequence.
        """
        unitig_graph_filename = os.path.join(os.path.dirname(__file__),
                                             'test_contig_placement_unitig_graph_4.gfa')
        unitig_graph = unicycler.string_graph.StringGraph(unitig_graph_filename)
        new_graph = unicycler.miniasm_assembly.place_contigs(self.working_dir, self.assembly_graph,
                                                             unitig_graph, 1, self.scoring_scheme,
                                                             self.seg_nums_to_bridge)
        self.assertTrue(len(new_graph.segments), 12)
        self.assertEqual(new_graph.get_total_segment_length(), 18391)

        merged_seqs = get_merged_string_graph_seqs(new_graph)
        self.assertTrue(sequences_match_some_rotation(merged_seqs[0],
                                                      unitig_graph.segments['1'].forward_sequence))
        self.assertTrue(sequences_match_some_rotation(merged_seqs[1],
                                                      unitig_graph.segments['2'].forward_sequence))
        self.assertTrue(sequences_match_some_rotation(merged_seqs[2],
                                                      unitig_graph.segments['3'].forward_sequence))
        self.assertTrue(sequences_match_some_rotation(merged_seqs[3],
                                                      unitig_graph.segments['4'].forward_sequence))
