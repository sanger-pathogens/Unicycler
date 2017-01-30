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
import unicycler.cpp_function_wrappers
import unicycler.read_ref
import unicycler.alignment


class TestSemiGlobalAlignment(unittest.TestCase):
    pass


class TestFullyGlobalAlignment(unittest.TestCase):
    pass


class TestPathAlignment(unittest.TestCase):
    pass


class TestMultipleSequenceAlignment(unittest.TestCase):

    def setUp(self):
        test_fastq = os.path.join(os.path.dirname(__file__), 'test_2.fastq')
        self.read_dict, self.read_names, _ = unicycler.read_ref.load_long_reads(test_fastq, 0)
        self.scoring_scheme = unicycler.alignment.AlignmentScoringScheme('3,-6,-5,-2')
        self.seqs = [self.read_dict[name].sequence for name in self.read_names]
        self.quals = [self.read_dict[name].qualities for name in self.read_names]
        self.original_seq = self.seqs[0]

    def test_consensus_with_subs(self):
        seqs = self.seqs[1:4]
        quals = self.quals[1:4]
        consensus, scores = unicycler.cpp_function_wrappers.consensus_alignment(seqs, quals,
                                                                                self.scoring_scheme)
        self.assertEqual(consensus, self.original_seq)

    def test_consensus_with_deletions(self):
        seqs = self.seqs[4:7]
        quals = self.quals[4:7]
        consensus, scores = unicycler.cpp_function_wrappers.consensus_alignment(seqs, quals,
                                                                                self.scoring_scheme)
        self.assertEqual(consensus, self.original_seq)

    def test_consensus_with_insertions(self):
        seqs = self.seqs[7:10]
        quals = self.quals[7:10]
        consensus, scores = unicycler.cpp_function_wrappers.consensus_alignment(seqs, quals,
                                                                                self.scoring_scheme)
        self.assertEqual(consensus, self.original_seq)

    def test_consensus_with_deletions_and_insertions(self):
        seqs = self.seqs[4:10]
        quals = self.quals[4:10]
        consensus, scores = unicycler.cpp_function_wrappers.consensus_alignment(seqs, quals,
                                                                                self.scoring_scheme)
        self.assertEqual(consensus, self.original_seq)

    def test_consensus_with_all(self):
        seqs = self.seqs[1:10]
        quals = self.quals[1:10]
        consensus, scores = unicycler.cpp_function_wrappers.consensus_alignment(seqs, quals,
                                                                                self.scoring_scheme)
        self.assertEqual(consensus, self.original_seq)

    def test_two_way_consensus(self):
        seqs = self.seqs[10:12]
        quals = self.quals[10:12]
        consensus, scores = unicycler.cpp_function_wrappers.consensus_alignment(seqs, quals,
                                                                                self.scoring_scheme)
        self.assertEqual(consensus, self.original_seq)

    def test_consensus_different_qualities(self):
        seqs = self.seqs[12:16]
        quals = self.quals[12:16]
        consensus, scores = unicycler.cpp_function_wrappers.consensus_alignment(seqs, quals,
                                                                                self.scoring_scheme)
        self.assertEqual(consensus, self.original_seq)
        self.assertEqual(scores[0], 1.0)
        self.assertTrue(scores[0] > scores[1])
        self.assertTrue(scores[1] > scores[2])
        self.assertTrue(scores[2] > scores[3])
