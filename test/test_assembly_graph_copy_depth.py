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
import unicycler.assembly_graph
import unicycler.assembly_graph_copy_depth


class TestCopyDepth(unittest.TestCase):

    def setUp(self):
        test_fastg = os.path.join(os.path.dirname(__file__), 'test_assembly_graph.fastg')
        self.graph = unicycler.assembly_graph.AssemblyGraph(test_fastg, 25, paths_file=None,
                                                            insert_size_mean=401,
                                                            insert_size_deviation=60)
        verbosity = 0
        unicycler.assembly_graph_copy_depth.determine_copy_depth(self.graph, verbosity)


    def test_single_copy_segments_1(self):
        """
        Tests large contigs from the largest replicon, all of which should have a copy depth of 1.
        """
        self.assertEqual(len(self.graph.copy_depths[141]), 1)
        self.assertEqual(len(self.graph.copy_depths[41]), 1)
        self.assertEqual(len(self.graph.copy_depths[306]), 1)
        self.assertEqual(len(self.graph.copy_depths[299]), 1)
        self.assertEqual(len(self.graph.copy_depths[125]), 1)
        self.assertEqual(len(self.graph.copy_depths[276]), 1)

    def test_single_copy_segments_2(self):
        """
        Tests single-copy contigs from a smaller, higher depth replicon. They should have a copy
        depth of 1, even though they are at a higher depth than the biggest replicon.
        """
        self.assertEqual(len(self.graph.copy_depths[272]), 1)
        self.assertEqual(len(self.graph.copy_depths[300]), 1)

    def test_merge_1(self):
        """
        Tests a particular part of the graph with merging and splitting.
        """
        self.assertEqual(len(self.graph.copy_depths[67]), 1)
        self.assertEqual(len(self.graph.copy_depths[165]), 1)
        self.assertEqual(len(self.graph.copy_depths[174]), 1)
        self.assertEqual(len(self.graph.copy_depths[242]), 1)
        self.assertEqual(len(self.graph.copy_depths[66]), 2)
        self.assertEqual(len(self.graph.copy_depths[65]), 3)
        self.assertEqual(len(self.graph.copy_depths[30]), 4)

    def test_merge_split_1(self):
        """
        Tests a particular part of the graph with merging and splitting.
        """
        self.assertEqual(len(self.graph.copy_depths[151]), 1)
        self.assertEqual(len(self.graph.copy_depths[114]), 1)
        self.assertEqual(len(self.graph.copy_depths[125]), 1)
        self.assertEqual(len(self.graph.copy_depths[152]), 2)
        self.assertEqual(len(self.graph.copy_depths[297]), 3)
        self.assertEqual(len(self.graph.copy_depths[55]), 1)
        self.assertEqual(len(self.graph.copy_depths[56]), 2)
        self.assertEqual(len(self.graph.copy_depths[222]), 3)
        self.assertEqual(len(self.graph.copy_depths[72]), 1)
        self.assertEqual(len(self.graph.copy_depths[137]), 2)
        self.assertEqual(len(self.graph.copy_depths[135]), 1)
        self.assertEqual(len(self.graph.copy_depths[136]), 1)

    def test_merge_split_2(self):
        """
        Tests a particular part of the graph with merging and splitting.
        """
        self.assertEqual(len(self.graph.copy_depths[271]), 1)
        self.assertEqual(len(self.graph.copy_depths[33]), 1)
        self.assertEqual(len(self.graph.copy_depths[232]), 2)
        self.assertEqual(len(self.graph.copy_depths[329]), 1)
        self.assertEqual(len(self.graph.copy_depths[330]), 1)
        self.assertEqual(len(self.graph.copy_depths[171]), 2)
        self.assertEqual(len(self.graph.copy_depths[172]), 1)
        self.assertEqual(len(self.graph.copy_depths[173]), 1)
        self.assertEqual(len(self.graph.copy_depths[309]), 2)
        self.assertEqual(len(self.graph.copy_depths[50]), 1)
        self.assertEqual(len(self.graph.copy_depths[308]), 3)
        self.assertEqual(len(self.graph.copy_depths[9]), 1)
        self.assertEqual(len(self.graph.copy_depths[10]), 2)
